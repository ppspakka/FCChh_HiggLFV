#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import argparse
from typing import Dict, List, Tuple, Optional

import ROOT

# ------------------------------
# Defaults and configuration
# ------------------------------

DEFAULT_IN_DIR = "./outdir"
DEFAULT_OUT_CARDS_DIR = "./datacards"
DEFAULT_OUT_ROOT = "./datacards/merged.root"
DEFAULT_CHANNEL = "bin1"

DEFAULT_INITIAL_HIST_NAME = "00_Initial_n_muons"
DEFAULT_FINAL_HIST_REGEX = r"^\d+_finalstate_nocut_m_collinear$"

# --- Configuration ---
JSON_PATH = "xsec_pb.json"
SIGNAL_TYPE = "ZH"  # Options: "ZH" or "VBF"

factors = {
    'ZToLL': 0.03363 * 2,    # Z->ee + Z->mumu
    'TauToLep': 0.1782 * 1,  # tau->e + tau->mu (consider opposite flavor only)
}

channel_map = {
    "etau": "ETauMu",
    "mutau": "MuTauE",
    "mue": "MuE"
}

def get_xsec(SIGNAL_TYPE):
    # --- Load JSON Data ---
    with open(JSON_PATH, "r") as jf:
        all_xsecs = json.load(jf)

    cross_sections_pb = {}
    signal_xsec_pb = {}

    # --- Parse JSON ---
    for key, xsec in all_xsecs.items():
        
        # Identify if the key is a signal (starts with ZH_ or VBF_)
        is_signal = key.startswith("ZH_") or key.startswith("VBF_")
        
        if is_signal:
            # Only process if it matches the desired SIGNAL_TYPE, ignore the other
            if key.startswith(SIGNAL_TYPE):
                # re.search finds the channel and mass regardless of the prefix (ZH_ll_ or VBF_)
                match = re.search(r"(etau|mutau)(\d+)", key)
                match_mutaue = re.search(r"mutau(\d+)", key)
                
                if match:
                    channel_raw, mass = match.groups()
                    channel_name = channel_map.get(channel_raw, channel_raw)
                    proc_name = f"H{channel_name}_LFV_{mass}"
                    
                    # Base scaling
                    scaled_xsec = xsec
                    
                    # Mass-dependent scaling
                    mass_int = int(mass)
                    if 110 <= mass_int <= 145:
                        scaled_xsec *= 1e-4
                    elif 145 < mass_int <= 190:
                        scaled_xsec *= 1e-2
                    elif 195 <= mass_int <= 240:
                        scaled_xsec *= 1.0
                        
                    if match_mutaue:
                        mue_proc_name = f"HMuE_LFV_{mass}"
                        signal_xsec_pb[mue_proc_name] = scaled_xsec

                    signal_xsec_pb[proc_name] = scaled_xsec * factors['TauToLep'] # Scale for tau->lepton branching
        else:
            # Treat all non-signal keys (zz_ll_tautau, zh_ll_ww, etc.) as backgrounds
            cross_sections_pb[key] = xsec

    # Combine backgrounds and the parsed/scaled signals
    cross_sections_pb.update(signal_xsec_pb)

    return cross_sections_pb


# Uncertainties framework (editable)
GLOBAL_UNC = {
    # applies to all processes (signal and background)
    "lumi": 1.02,
}

BACKGROUND_UNC_NAME = "bkg_unc"
BACKGROUND_UNC_VALUE = 1.3 # 30% uncertainty on backgrounds

UNCERTAINTIES = {
    "lumi": 1.02, # This indicate apply globally
    "example_proc_unc": {
        "proc1": 1.05,
        "proc2": 1.10,
    },
    "DY" : {
        "DY0j": 1.3,
        "DY1j": 1.3,
    },
    "ttbar" : {
        "ttbar": 1.3,
    },
    "tW" : {
        "tW": 1.3,
        "tbarW": 1.3,
    },
    "SM_ggH" : {
        "SM_ggH_tautau": 1.3,
        "SM_ggH_WW": 1.3,
    },
    "SM_VBF" : {
        "SM_VBFH_tautau": 1.3,
        "SM_VBFH_WW": 1.3,
    },
    "WW" : {
        "WW": 1.3,
    },
}

# ------------------------------
# Helpers
# ------------------------------

def die(msg: str, code: int = 2):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)

def info(msg: str):
    print(f"[INFO] {msg}")

def is_root_file(fname: str) -> bool:
    return fname.endswith(".root") and os.path.isfile(fname)

def classify_file(fname: str) -> Optional[Tuple[str, str]]:
    """
    Return (kind, process_name)
    kind: "signal" or "background"
    process_name: stripped name (remove prefix signal_/background_ and .root)
    """
    base = os.path.basename(fname)
    if base.startswith("signal_") and base.endswith(".root"):
        return "signal", base[len("signal_"):-len(".root")]
    if base.startswith("background_") and base.endswith(".root"):
        return "background", base[len("background_"):-len(".root")]
    return None

def get_all_objects(tdir: ROOT.TDirectory):
    """Yield (path, obj) for all objects recursively under tdir."""
    for key in tdir.GetListOfKeys():
        obj = key.ReadObj()
        name = obj.GetName()
        path = f"{tdir.GetPath()}/{name}"
        if obj.InheritsFrom("TDirectory"):
            yield from get_all_objects(obj)
        else:
            yield path, obj

def find_hist_by_name(tfile: ROOT.TFile, target_name: str) -> Optional[ROOT.TH1]:
    """Find first TH1 with exact name target_name anywhere in file."""
    for _, obj in get_all_objects(tfile):
        if obj.InheritsFrom("TH1") and obj.GetName() == target_name:
            return obj
    return None

def find_unique_final_hist(tfile: ROOT.TFile, pattern: re.Pattern) -> ROOT.TH1:
    matches = []
    for _, obj in get_all_objects(tfile):
        if obj.InheritsFrom("TH1"):
            name = obj.GetName()
            if pattern.match(name):
                matches.append(obj)
    if len(matches) == 0:
        die(f"No histogram matching final pattern in file: {tfile.GetName()}")
    if len(matches) > 1:
        die(f"Multiple histograms match final pattern in file: {tfile.GetName()}")
    return matches[0]

def compute_weight(xsec_pb: float, lumi_pb: float, total_events: float) -> float:
    if total_events <= 0:
        die("Total events in initial histogram is zero; cannot compute weight.")
    return (xsec_pb * lumi_pb) / total_events

def clone_and_scale_hist(h: ROOT.TH1, new_name: str, scale: float) -> ROOT.TH1:
    out = h.Clone(new_name)
    out.SetDirectory(0)  # detach from input file
    out.Scale(scale)
    return out

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def write_merged_root(
    out_root_path: str,
    channel: str,
    proc_to_hist: Dict[str, ROOT.TH1],
    data_obs_hist: Optional[ROOT.TH1] = None
):
    ensure_dir(os.path.dirname(out_root_path) or ".")
    f = ROOT.TFile(out_root_path, "RECREATE")
    if not f or f.IsZombie():
        die(f"Cannot create output ROOT: {out_root_path}")

    # Create per-channel directory
    ch_dir = f.mkdir(channel)
    ch_dir.cd()

    # Write one hist per process with the process name
    for proc, hist in proc_to_hist.items():
        hwrite = hist.Clone(proc)
        hwrite.SetDirectory(ch_dir)
        
        # for ibin in range(1, hwrite.GetNbinsX() + 1):
        #     hwrite.SetBinError(ibin, 0.0)
        
        
        hwrite.Write()

    # Write data_obs if provided (sum of all backgrounds)
    if data_obs_hist is not None:
        dobj = data_obs_hist.Clone("data_obs")
        dobj.SetDirectory(ch_dir)
        
        # for ibin in range(1, dobj.GetNbinsX() + 1):
        #     dobj.SetBinError(ibin, 0.0)
        
        dobj.Write()

    f.Write()
    f.Close()
    n_total = len(proc_to_hist) + (1 if data_obs_hist is not None else 0)
    info(f"Wrote merged ROOT with {n_total} histograms at {out_root_path}")

def maybe_rebin_variable(hist: ROOT.TH1, bin_edges: Optional[List[float]]) -> ROOT.TH1:
    """
    Optional variable rebinning support (not applied by default).
    Pass bin_edges as a sorted list of edges if rebin is desired.
    """
    if not bin_edges:
        return hist
    import array
    edges = array.array('d', bin_edges)
    rebinned = hist.Rebin(len(edges) - 1, hist.GetName() + "_rebin", edges)
    return rebinned

def get_yield_around_mass(hist: ROOT.TH1, mass: float, window: float) -> float:
    """
    For counting experiment
    sum bin contents in [mass - window, mass + window]
    """
    yield_sum = 0.0
    # for ibin in range(1, hist.GetNbinsX() + 1):
    #     bin_center = hist.GetBinCenter(ibin)
    #     if mass - window <= bin_center <= mass + window:
    #         yield_sum += hist.GetBinContent(ibin)
    # No attribute GetNbinsX() error fix
    nbins = hist.GetXaxis().GetNbins()
    for ibin in range(1, nbins + 1):
        bin_center = hist.GetBinCenter(ibin)
        if mass - window <= bin_center <= mass + window:
            yield_sum += hist.GetBinContent(ibin)
    return yield_sum

def write_unc_line(procs: List[str]): # procs for the indexing
    global UNCERTAINTIES
    lines = []
    for unc_name, unc_val in UNCERTAINTIES.items():
        # Global uncertainty
        if isinstance(unc_val, float):
            parts = [unc_name, "lnN"]
            for _ in procs:
                parts.append(f"{unc_val:.3g}")
            lines.append(parts)
        # Process-specific uncertainties
        elif isinstance(unc_val, dict):
            parts = [unc_name, "lnN"]
            for p in procs:
                if p in unc_val:
                    parts.append(f"{unc_val[p]:.3g}")
                else:
                    parts.append("-")
            lines.append(parts)
    return lines

def write_datacard_for_signal(
    signal_proc: str,
    bkg_procs: List[str],
    channel: str,
    shapes_root_path: str,
    out_dir: str,
    global_unc: Dict[str, float],
    bkg_unc_name: str,
    bkg_unc_value: float,
    observation: int = 0,
    mode: str = "counting", # shape or counting
):
    ensure_dir(out_dir)

    # Order: signal first, then backgrounds (sorted for reproducibility)
    procs = [signal_proc] + sorted(bkg_procs)
    nprocs = len(procs)

    # Process ID line: 0 for signal, 1.. for bkgs
    proc_ids = [0] + list(range(1, nprocs))

    # Build uncertainty lines
    # GLOBAL_UNC applies to all processes
    global_unc_lines = []
    for unc_name, unc_val in global_unc.items():
        parts = [unc_name, "lnN"]
        for _ in procs:
            parts.append(f"{unc_val:.3g}")
        global_unc_lines.append(parts)

    # Background uncertainty applies only to backgrounds
    bkg_unc_line = [bkg_unc_name, "lnN"]
    for p in procs:
        if p == signal_proc:
            bkg_unc_line.append("-")
        else:
            bkg_unc_line.append(f"{bkg_unc_value:.3g}")

    # Formatting helpers
    def fmt_row(cells: List[str], widths=None) -> str:
        if widths is None:
            return "  ".join(str(c) for c in cells)
        return "  ".join(str(c).ljust(w) for c, w in zip(cells, widths))

    # Compute column widths for alignment
    col_headers = ["bin"] + procs
    widths = [max(len(str(x)), 8) for x in col_headers]

    card_path = os.path.join(out_dir, f"datacard_{signal_proc}.txt")
    
    # TEMP FIX: for shapes_root_path, extract the last file name 
    file_path = shapes_root_path.split('/')[-1]
    shapes_root_path = f"{file_path}"
    
    if mode == "shape":
        with open(card_path, "w") as dc:
            dc.write("imax 1 number of channels\n")
            dc.write("jmax * number of backgrounds\n")
            dc.write("kmax * number of nuisance parameters (sources of systematical uncertainties)\n")
            dc.write("------------\n")
            dc.write(f"shapes data_obs {channel} {shapes_root_path} {channel}/data_obs\n")
            dc.write(f"shapes *        {channel} {shapes_root_path} {channel}/$PROCESS\n")
            dc.write("------------\n")
            dc.write(f"bin {channel}\n")
            dc.write(f"observation {observation}\n")
            dc.write("------------\n")

            # bin line for each process
            dc.write(fmt_row(["bin"] + [channel] * nprocs, widths) + "\n")
            # process names
            dc.write(fmt_row(["process"] + procs, widths) + "\n")
            # process ids
            dc.write(fmt_row(["process"] + [str(i) for i in proc_ids], widths) + "\n")
            # rates from shapes
            dc.write(fmt_row(["rate"] + ["-1"] * nprocs, widths) + "\n")
            dc.write("------------\n")

            # Uncertainties
            # for parts in global_unc_lines:
            #     dc.write(fmt_row(parts) + "\n")
            # dc.write(fmt_row(bkg_unc_line) + "\n")
            
            # Use new uncertainty writing function
            unc_lines = write_unc_line(procs)
            for parts in unc_lines:
                dc.write(fmt_row(parts) + "\n")
            
            dc.write("\n")

            # autoMCStats and rateParam/freeze
            # dc.write("autoMCStats 1\n")
            dc.write("lumiscale rateParam     *           *           1.0\n")
            dc.write("nuisance edit freeze lumiscale\n")
    elif mode == "counting":
        realfile_path = os.path.join(out_dir, file_path)
        
        signal_hist = None
        file = ROOT.TFile.Open(realfile_path)
        bin1_dir = file.Get(channel)
        for key in bin1_dir.GetListOfKeys():
            if key.GetName() == signal_proc:
                signal_hist = bin1_dir.Get(key.GetName())
            
        signal_yield = get_yield_around_mass(
            signal_hist,
            mass=float(signal_proc.split('_')[-1]),
            window=10.0,
        )
        
        background_hist = {}
        for key in bin1_dir.GetListOfKeys():
            bkg_name = key.GetName()
            if bkg_name in bkg_procs:
                background_hist[bkg_name] = bin1_dir.Get(bkg_name)
        bkg_yields = []
        for bkg in bkg_procs:
            byield = get_yield_around_mass(
                background_hist[bkg],
                mass=float(signal_proc.split('_')[-1]),
                window=10.0,
            )
            if byield == 0.0:
                warn(f"Background '{bkg}' has zero yield around mass for datacard {signal_proc} set to 1e-12 to avoid issues.")
                byield = 1e-12
            bkg_yields.append(byield)
        # Print all yields
        info(f"Yields for datacard {signal_proc}: signal={signal_yield:.6g}, backgrounds={[f'{by:.6g}' for by in bkg_yields]}")
        
        # Calculate the observation as sum of backgrounds (rounded)
        if observation == 0:
            observation = int(round(sum(bkg_yields)))
            info(f"Setting observation to sum of backgrounds: {observation}")
        
        with open(card_path, "w") as dc:
            dc.write("imax 1 number of channels\n")
            dc.write("jmax * number of backgrounds\n")
            dc.write("kmax * number of nuisance parameters (sources of systematical uncertainties)\n")
            dc.write("------------\n")
            dc.write(f"bin {channel}\n")
            dc.write(f"observation {observation}\n")
            dc.write("------------\n")

            # bin line for each process
            dc.write(fmt_row(["bin"] + [channel] * nprocs, widths) + "\n")
            # process names
            dc.write(fmt_row(["process"] + procs, widths) + "\n")
            # process ids
            dc.write(fmt_row(["process"] + [str(i) for i in proc_ids], widths) + "\n")
            # rates from counting
            rate_strs = [f"{signal_yield:.6g}"] + [f"{by:.6g}" for by in bkg_yields]
            dc.write(fmt_row(["rate"] + rate_strs, widths) + "\n")
            dc.write("------------\n")

            # Uncertainties
            # Use new uncertainty writing function
            unc_lines = write_unc_line(procs)
            for parts in unc_lines:
                dc.write(fmt_row(parts) + "\n")
            dc.write("\n")

            # autoMCStats and rateParam/freeze
            # dc.write("autoMCStats 1\n")
            dc.write("lumiscale rateParam     *           *           1.0\n")
            dc.write("nuisance edit freeze lumiscale\n")
    info(f"Wrote datacard: {card_path}")

# ------------------------------
# Main workflow
# ------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate shapes and Combine datacards for FCCee HLFV study.")
    parser.add_argument("--in-dir", default=DEFAULT_IN_DIR, help="Input directory containing ROOT files (signal_*.root, background_*.root)")
    parser.add_argument("--out-root", default=DEFAULT_OUT_ROOT, help="Output merged shapes ROOT path")
    parser.add_argument("--out-cards", default=DEFAULT_OUT_CARDS_DIR, help="Output directory for datacards")
    parser.add_argument("--channel", default=DEFAULT_CHANNEL, help="Channel name (directory inside shapes ROOT)")
    parser.add_argument("--lumi-pb", type=float, required=True, help="Target integrated luminosity in pb^-1 (e.g., 1e7 for 10 ab^-1)")
    parser.add_argument("--initial-hist", default=DEFAULT_INITIAL_HIST_NAME, help="Name of the initial histogram for total entry counting")
    parser.add_argument("--final-hist-pattern", default=DEFAULT_FINAL_HIST_REGEX, help="Regex for selecting the final histogram")
    parser.add_argument("--xsec-json", default=None, help="Optional JSON file to override/extend cross sections dict")
    parser.add_argument("--rebin-json", default=None, help="Optional JSON with variable bin edges: {'bin1': [edges...]}")
    parser.add_argument("--signal-type", type=str, choices=["ZH", "VBF"], default=SIGNAL_TYPE, help="Type of signal to process (ZH or VBF)")
    args = parser.parse_args()

    in_dir = args.in_dir
    out_root = args.out_root
    out_cards = args.out_cards
    channel = args.channel
    lumi_pb = args.lumi_pb

    if not os.path.isdir(in_dir):
        die(f"Input directory not found: {in_dir}")

    # Cross-section map
    cross_sections_pb = get_xsec(args.signal_type)
    xsec_map = dict(cross_sections_pb)
    if args.xsec_json:
        if not os.path.isfile(args.xsec_json):
            die(f"xsec-json not found: {args.xsec_json}")
        with open(args.xsec_json) as jf:
            override = json.load(jf)
        xsec_map.update(override)
        info(f"Loaded cross-section overrides for {len(override)} processes")

    # Rebin config (not applied by default, just prepared)
    rebin_edges = None
    if args.rebin_json:
        if not os.path.isfile(args.rebin_json):
            die(f"rebin-json not found: {args.rebin_json}")
        with open(args.rebin_json) as jf:
            rebin_cfg = json.load(jf)
        rebin_edges = rebin_cfg.get(channel, None)
        info(f"Loaded rebin config for channel {channel}: {('enabled' if rebin_edges else 'none')}")

    # Collect files
    files = [
        os.path.join(in_dir, f)
        for f in os.listdir(in_dir)
        if is_root_file(os.path.join(in_dir, f))
    ]
    if not files:
        die(f"No ROOT files found in {in_dir}")

    signals: List[str] = []
    backgrounds: List[str] = []
    proc_to_hist: Dict[str, ROOT.TH1] = {}

    final_pattern = re.compile(args.final_hist_pattern)

    for fpath in sorted(files):
        kind_proc = classify_file(fpath)
        if kind_proc is None:
            warn(f"Skipping non-matching filename (expect signal_*.root or background_*.root): {os.path.basename(fpath)}")
            continue

        kind, proc = kind_proc
        # Ensure cross-section exists
        if proc not in xsec_map:
            die(f"Missing cross-section for process '{proc}'. Add to cross_sections_pb or provide --xsec-json.")

        tf = ROOT.TFile.Open(fpath)
        if not tf or tf.IsZombie():
            die(f"Failed to open ROOT file: {fpath}")

        # Find initial and final hists
        initial_hist = find_hist_by_name(tf, args.initial_hist)
        if initial_hist is None:
            tf.Close()
            die(f"Initial histogram '{args.initial_hist}' not found in file: {fpath}")

        final_hist = find_unique_final_hist(tf, final_pattern)

        total_events = float(initial_hist.GetEntries())
        # per user spec, entries is the right total (raw counting), not Integral()
        if total_events <= 0:
            tf.Close()
            die(f"Initial histogram has zero entries in file: {fpath}")

        xsec = float(xsec_map[proc])  # pb
        weight = compute_weight(xsec, lumi_pb, total_events)

        # Efficiency and yield for logging
        passed = float(final_hist.GetEntries())
        eff = passed / total_events if total_events > 0 else 0.0
        # Weight applies per-event to final hist integrals; yield ~ integral * weight if integral counts entries
        est_yield = passed * weight

        # Prepare hist for merged root
        h = clone_and_scale_hist(final_hist, proc, weight)
        if rebin_edges:
            h = maybe_rebin_variable(h, rebin_edges)

        tf.Close()

        proc_to_hist[proc] = h
        if kind == "signal":
            signals.append(proc)
        else:
            backgrounds.append(proc)

        info(f"Process: {proc:20s} kind={kind:10s} xsec={xsec:.6g} pb  total={total_events:.0f}  eff={eff:.4f}  weight={weight:.6g}  yield~{est_yield:.6g}")

    if not signals:
        die("No signal_*.root found.")
    if not backgrounds:
        die("No background_*.root found.")

    # Build data_obs as the sum of all background histograms
    data_obs = None
    for i, b in enumerate(backgrounds):
        hb = proc_to_hist[b]
        if i == 0:
            data_obs = hb.Clone("data_obs")
            data_obs.SetDirectory(0)
        else:
            data_obs.Add(hb)

    # Write merged ROOT shapes (including data_obs)
    write_merged_root(out_root, channel, proc_to_hist, data_obs_hist=data_obs)

    # Build PROC_UNC dynamically for backgrounds (for future extension if needed)
    # Currently used only to materialize the "bkg_unc" line.
    proc_unc = {b: {BACKGROUND_UNC_NAME: BACKGROUND_UNC_VALUE} for b in backgrounds}

    # Generate one datacard per signal
    for sig in sorted(signals):
        write_datacard_for_signal(
            signal_proc=sig,
            bkg_procs=backgrounds,
            channel=channel,
            shapes_root_path=out_root,
            out_dir=out_cards,
            global_unc=GLOBAL_UNC,
            bkg_unc_name=BACKGROUND_UNC_NAME,
            bkg_unc_value=BACKGROUND_UNC_VALUE,
            observation=0,  # projection study
        )

    info(f"Done. Datacards written to {out_cards}")

if __name__ == "__main__":
    main()
