#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import argparse
import array
from typing import Dict, List, Tuple, Optional

import ROOT

# ------------------------------
# Defaults and configuration
# ------------------------------

DEFAULT_IN_DIR = "./outdir"
DEFAULT_OUT_CARDS_DIR = "./datacards"
DEFAULT_CHANNEL = "bin1"

DEFAULT_INITIAL_HIST_NAME = "00_Initial_n_muons"
DEFAULT_FINAL_HIST_REGEX = r"^\d+_finalstate_nocut_m_collinear$"

# --- Configuration ---
JSON_PATH = "xsec_pb.json"
SIGNAL_TYPE = "ggH"  # Options: "ZH" or "VBF"
WINDOW = 50.0

# Unequal binning definition for Shape Analysis (edges in GeV)
# The original contains 1400 bins from 0 to 1400 GeV.
# This variable array merges bins (narrow around peak, wider at tails).
UNEQUAL_BIN_EDGES = [
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
    110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
    250.0, 300.0, 400.0, 500.0, 700.0, 1000.0, 1400.0
]
lowmass_bin_edges = list(range(75, 926, 50)) + [1000, 1100, 1200, 1400]
highmass_bin_edges = [0, 75, 125, 175, 225, 275, 325, 375, 525, 675, 825, 975, 1400]

factors = {
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
        is_signal = key.startswith("ZH_") or key.startswith("VBF_") or key.startswith("ggH_")
        
        if is_signal:
            # Only process if it matches the desired SIGNAL_TYPE, ignore the other
            if key.startswith(SIGNAL_TYPE):
                signal_xsec_pb[key] = xsec * factors['TauToLep']
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
    base = os.path.basename(fname)
    if base.startswith("signal_") and base.endswith(".root"):
        return "signal", base[len("signal_"):-len(".root")]
    if base.startswith("background_") and base.endswith(".root"):
        return "background", base[len("background_"):-len(".root")]
    return None

def get_all_objects(tdir: ROOT.TDirectory):
    for key in tdir.GetListOfKeys():
        obj = key.ReadObj()
        name = obj.GetName()
        path = f"{tdir.GetPath()}/{name}"
        if obj.InheritsFrom("TDirectory"):
            yield from get_all_objects(obj)
        else:
            yield path, obj

def find_hist_by_name(tfile: ROOT.TFile, target_name: str) -> Optional[ROOT.TH1]:
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

def write_shapes_root(
    out_root_path: str,
    channel: str,
    signal_proc: str,
    signal_hist: ROOT.TH1,
    bkg_hists_dict: Dict[str, ROOT.TH1],
    data_obs_hist: Optional[ROOT.TH1] = None
):
    """Write a specific root file containing the signal and all backgrounds for shape analysis."""
    ensure_dir(os.path.dirname(out_root_path) or ".")
    f = ROOT.TFile(out_root_path, "RECREATE")
    if not f or f.IsZombie():
        die(f"Cannot create output ROOT: {out_root_path}")

    ch_dir = f.mkdir(channel)
    ch_dir.cd()

    # Write signal
    h_sig = signal_hist.Clone(signal_proc)
    h_sig.SetDirectory(ch_dir)
    h_sig.Write()

    # Write backgrounds
    for bkg_name, h_bkg in bkg_hists_dict.items():
        h_b = h_bkg.Clone(bkg_name)
        h_b.SetDirectory(ch_dir)
        h_b.Write()

    # Write data_obs
    if data_obs_hist is not None:
        dobj = data_obs_hist.Clone("data_obs")
        dobj.SetDirectory(ch_dir)
        dobj.Write()

    f.Write()
    f.Close()
    info(f"Wrote shapes ROOT for {signal_proc} at {out_root_path}")

def maybe_rebin_variable(hist: ROOT.TH1, bin_edges: Optional[List[float]]) -> ROOT.TH1:
    """Optional variable rebinning support."""
    if not bin_edges:
        return hist
    edges = array.array('d', bin_edges)
    rebinned = hist.Rebin(len(edges) - 1, hist.GetName() + "_rebin", edges)
    rebinned.SetDirectory(0)
    return rebinned

def get_yield_around_mass(hist: ROOT.TH1, mass: float, window: float) -> float:
    yield_sum = 0.0
    nbins = hist.GetXaxis().GetNbins()
    for ibin in range(1, nbins + 1):
        bin_center = hist.GetBinCenter(ibin)
        if mass - window <= bin_center <= mass + window:
            yield_sum += hist.GetBinContent(ibin)
    return yield_sum

def write_unc_line(procs: List[str]):
    global UNCERTAINTIES
    lines = []
    for unc_name, unc_val in UNCERTAINTIES.items():
        if isinstance(unc_val, float):
            parts = [unc_name, "lnN"]
            for _ in procs:
                parts.append(f"{unc_val:.3g}")
            lines.append(parts)
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
    shapes_root_name: str,
    out_dir: str,
    global_unc: Dict[str, float],
    bkg_unc_name: str,
    bkg_unc_value: float,
    mode: str,
    signal_hist: ROOT.TH1,
    bkg_hists: Dict[str, ROOT.TH1],
    observation: int = 0,
):
    ensure_dir(out_dir)

    procs = [signal_proc] + sorted(bkg_procs)
    nprocs = len(procs)
    proc_ids = [0] + list(range(1, nprocs))

    def fmt_row(cells: List[str], widths=None) -> str:
        if widths is None:
            return "  ".join(str(c) for c in cells)
        return "  ".join(str(c).ljust(w) for c, w in zip(cells, widths))

    col_headers = ["bin"] + procs
    widths = [max(len(str(x)), 8) for x in col_headers]

    card_path = os.path.join(out_dir, f"datacard_{signal_proc}.txt")
    
    if mode == "shape":
        with open(card_path, "w") as dc:
            dc.write("imax 1 number of channels\n")
            dc.write("jmax * number of backgrounds\n")
            dc.write("kmax * number of nuisance parameters (sources of systematical uncertainties)\n")
            dc.write("------------\n")
            dc.write(f"shapes data_obs {channel} {shapes_root_name} {channel}/data_obs\n")
            dc.write(f"shapes *        {channel} {shapes_root_name} {channel}/$PROCESS\n")
            dc.write("------------\n")
            dc.write(f"bin {channel}\n")
            dc.write(f"observation {observation}\n")
            dc.write("------------\n")
            dc.write(fmt_row(["bin"] + [channel] * nprocs, widths) + "\n")
            dc.write(fmt_row(["process"] + procs, widths) + "\n")
            dc.write(fmt_row(["process"] + [str(i) for i in proc_ids], widths) + "\n")
            dc.write(fmt_row(["rate"] + ["-1"] * nprocs, widths) + "\n")
            dc.write("------------\n")

            unc_lines = write_unc_line(procs)
            for parts in unc_lines:
                dc.write(fmt_row(parts) + "\n")
            
            dc.write("\n")
            dc.write("lumiscale rateParam     *           *           1.0\n")
            dc.write("nuisance edit freeze lumiscale\n")
            
    elif mode == "counting":
        signal_yield = get_yield_around_mass(
            signal_hist,
            mass=float(signal_proc.split('_')[-1]),
            window=WINDOW,
        )
        
        bkg_yields = []
        for bkg in bkg_procs:
            byield = get_yield_around_mass(
                bkg_hists[bkg],
                mass=float(signal_proc.split('_')[-1]),
                window=WINDOW,
            )
            if byield == 0.0:
                warn(f"Background '{bkg}' has zero yield around mass for datacard {signal_proc} set to 1e-12 to avoid issues.")
                byield = 1e-12
            bkg_yields.append(byield)
            
        info(f"Yields for datacard {signal_proc}: signal={signal_yield:.6g}, backgrounds={[f'{by:.6g}' for by in bkg_yields]}")
        
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
            dc.write(fmt_row(["bin"] + [channel] * nprocs, widths) + "\n")
            dc.write(fmt_row(["process"] + procs, widths) + "\n")
            dc.write(fmt_row(["process"] + [str(i) for i in proc_ids], widths) + "\n")
            
            rate_strs = [f"{signal_yield:.6g}"] + [f"{by:.6g}" for by in bkg_yields]
            dc.write(fmt_row(["rate"] + rate_strs, widths) + "\n")
            dc.write("------------\n")

            unc_lines = write_unc_line(procs)
            for parts in unc_lines:
                dc.write(fmt_row(parts) + "\n")
            dc.write("\n")
            dc.write("lumiscale rateParam     *           *           1.0\n")
            dc.write("nuisance edit freeze lumiscale\n")
            
    info(f"Wrote datacard: {card_path}")

# ------------------------------
# Main workflow
# ------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate shapes and Combine datacards for FCCee HLFV study.")
    parser.add_argument("--mode", choices=["shape", "counting"], default="shape", help="Analysis mode: generate shapes or use counting experiment")
    parser.add_argument("--in-dir", default=DEFAULT_IN_DIR, help="Input directory containing ROOT files (signal_*.root, background_*.root)")
    parser.add_argument("--out-cards", default=DEFAULT_OUT_CARDS_DIR, help="Output directory for datacards and per-signal roots")
    parser.add_argument("--channel", default=DEFAULT_CHANNEL, help="Channel name (directory inside shapes ROOT)")
    parser.add_argument("--lumi-pb", type=float, required=True, help="Target integrated luminosity in pb^-1 (e.g., 1e7 for 10 ab^-1)")
    parser.add_argument("--initial-hist", default=DEFAULT_INITIAL_HIST_NAME, help="Name of the initial histogram for total entry counting")
    parser.add_argument("--final-hist-pattern", default=DEFAULT_FINAL_HIST_REGEX, help="Regex for selecting the final histogram")
    parser.add_argument("--xsec-json", default=None, help="Optional JSON file to override/extend cross sections dict")
    parser.add_argument("--rebin-json", default=None, help="Optional JSON with variable bin edges: {'bin1': [edges...]}")
    parser.add_argument("--signal-type", type=str, choices=["ZH", "VBF"], default=SIGNAL_TYPE, help="Type of signal to process (ZH or VBF)")
    parser.add_argument("--out-root", default=None, help="Old variable (unused), if input is provided, it will be ignored and a warning will be issued.")
    parser.add_argument("--bin-type", choices=["lowmass", "highmass"], default="lowmass", help="Predefined binning scheme to use instead of custom edges (lowmass or highmass)")
    args = parser.parse_args()

    in_dir = args.in_dir
    out_cards = args.out_cards
    channel = args.channel
    lumi_pb = args.lumi_pb
    
    if args.out_root is not None:
        warn("--out-root is no longer used and will be ignored. Shapes ROOT files will be written to the same directory as datacards with naming datacard_<signal>.root.")

    if not os.path.isdir(in_dir):
        die(f"Input directory not found: {in_dir}")

    cross_sections_pb = get_xsec(args.signal_type)
    xsec_map = dict(cross_sections_pb)
    if args.xsec_json:
        if not os.path.isfile(args.xsec_json):
            die(f"xsec-json not found: {args.xsec_json}")
        with open(args.xsec_json) as jf:
            override = json.load(jf)
        xsec_map.update(override)
        info(f"Loaded cross-section overrides for {len(override)} processes")

    # Determine rebin edges (JSON override > Default unequal)
    rebin_edges = None
    if args.rebin_json:
        if not os.path.isfile(args.rebin_json):
            die(f"rebin-json not found: {args.rebin_json}")
        with open(args.rebin_json) as jf:
            rebin_cfg = json.load(jf)
        rebin_edges = rebin_cfg.get(channel, None)
    elif args.bin_type == "lowmass":
        rebin_edges = lowmass_bin_edges
    elif args.bin_type == "highmass":
        rebin_edges = highmass_bin_edges

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
        if proc not in xsec_map:
            die(f"Missing cross-section for process '{proc}'. Add to cross_sections_pb or provide --xsec-json.")

        tf = ROOT.TFile.Open(fpath)
        if not tf or tf.IsZombie():
            die(f"Failed to open ROOT file: {fpath}")

        initial_hist = find_hist_by_name(tf, args.initial_hist)
        if initial_hist is None:
            tf.Close()
            die(f"Initial histogram '{args.initial_hist}' not found in file: {fpath}")

        final_hist = find_unique_final_hist(tf, final_pattern)

        total_events = float(initial_hist.GetEntries())
        if total_events <= 0:
            tf.Close()
            die(f"Initial histogram has zero entries in file: {fpath}")

        xsec = float(xsec_map[proc])  # pb
        weight = compute_weight(xsec, lumi_pb, total_events)

        passed = float(final_hist.GetEntries())
        eff = passed / total_events if total_events > 0 else 0.0
        est_yield = passed * weight

        h = clone_and_scale_hist(final_hist, proc, weight)
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

    # Build overarching unbinned data_obs
    data_obs = None
    for i, b in enumerate(backgrounds):
        hb = proc_to_hist[b]
        if i == 0:
            data_obs = hb.Clone("data_obs")
            data_obs.SetDirectory(0)
        else:
            data_obs.Add(hb)
    # Generate one datacard and root file per signal
    for sig in sorted(signals):
        sig_hist = proc_to_hist[sig]
        bkg_dict = {b: proc_to_hist[b] for b in backgrounds}

        # Apply Rebinning ONLY for shape analysis
        if args.mode == "shape" and rebin_edges:
            sig_h_use = maybe_rebin_variable(sig_hist, rebin_edges)
            bkg_dict_use = {b: maybe_rebin_variable(h, rebin_edges) for b, h in bkg_dict.items()}
            data_obs_use = maybe_rebin_variable(data_obs, rebin_edges)
        else:
            sig_h_use = sig_hist
            bkg_dict_use = bkg_dict
            data_obs_use = data_obs

        # === NEW: Filter out empty histograms to prevent "Null norm" error ===
        
        # 1. Check if the signal is empty. If it is, Combine cannot calculate a limit.
        if sig_h_use.Integral() <= 0:
            warn(f"Signal '{sig}' has zero total yield. Skipping datacard generation for this mass.")
            continue

        # 2. Filter out any backgrounds with exactly 0 yield
        active_backgrounds = []
        active_bkg_dict = {}
        for b_name, h_bkg in bkg_dict_use.items():
            if h_bkg.Integral() > 0:
                active_backgrounds.append(b_name)
                active_bkg_dict[b_name] = h_bkg
            else:
                warn(f"Background '{b_name}' has zero yield. Excluding from {sig} datacard.")

        # Set up per-signal naming
        root_filename = f"datacard_{sig}.root"
        out_root_path = os.path.join(out_cards, root_filename)

        # Write shapes ROOT file using ONLY the non-empty backgrounds
        write_shapes_root(out_root_path, channel, sig, sig_h_use, active_bkg_dict, data_obs_use)

        write_datacard_for_signal(
            signal_proc=sig,
            bkg_procs=active_backgrounds,  # <--- Pass the filtered list
            channel=channel,
            shapes_root_name=root_filename,
            out_dir=out_cards,
            global_unc=GLOBAL_UNC,
            bkg_unc_name=BACKGROUND_UNC_NAME,
            bkg_unc_value=BACKGROUND_UNC_VALUE,
            mode=args.mode,
            signal_hist=sig_h_use,
            bkg_hists=active_bkg_dict,     # <--- Pass the filtered dict
            observation=-1, 
        )

    info(f"Done. Datacards written to {out_cards}")

if __name__ == "__main__":
    main()