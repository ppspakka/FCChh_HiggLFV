import os
import re
import json
import subprocess
import argparse
import concurrent.futures
import logging
import shutil
from pathlib import Path

# Defaults
DEFAULT_LUMI = "1"
DEFAULT_SEED = "1"
DEFAULT_OUTPUT_JSON = "limits.json"
DEFAULT_THREADS = 4

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def get_prefix():
    status_dir = Path("status")
    queue = len(list(status_dir.glob("*.Queue")))
    running = len(list(status_dir.glob("*.Running")))
    complete = len(list(status_dir.glob("*.Complete")))
    return f"[Queue {queue}][Running {running}][Complete {complete}] "

def estimate_rmax(datacard: Path, safety_factor: float = 10.0) -> float:
    """
    Parses the datacard to estimate rMax based on stat uncertainty + a safety factor
    for future systematic uncertainties.
    """
    try:
        with open(datacard, 'r') as f:
            lines = f.readlines()
        
        processes = []
        rates = []
        
        for line in lines:
            parts = line.strip().split()
            if not parts: continue
            
            # Identify process IDs (e.g. 0 for signal, 1, 2... for bkg)
            if parts[0] == "process":
                if all(p.lstrip('-').isdigit() for p in parts[1:]):
                    processes = [int(p) for p in parts[1:]]
            # Identify expected yields
            elif parts[0] == "rate":
                rates = [float(p) for p in parts[1:]]
                
        if not processes or not rates:
            return 0.1 # Fallback 
            
        signal_rate = sum(r for p, r in zip(processes, rates) if p <= 0)
        bkg_rate = sum(r for p, r in zip(processes, rates) if p > 0)
        
        if signal_rate <= 0 or bkg_rate <= 0: 
            return 0.1 # Fallback
        
        # 95% CL expected stat limit ~ 1.96 * sqrt(B) / S
        r_expected = 1.96 * (bkg_rate ** 0.5) / signal_rate
        r_max = r_expected * safety_factor
        
        # Keep it to a sensible precision
        return float(f"{r_max:.4g}")
    except Exception as e:
        logger.warning(f"Failed to estimate rMax for {datacard.name}, using default 0.1. Error: {e}")
        return 0.1

def parse_limits(log: str):
    """
    Extracts all quantiles and observed limit from the AsymptoticLimits stdout.
    """
    results = {}
    for line in log.splitlines():
        if "Observed Limit:" in line:
            m = re.search(r"r\s*<\s*([0-9\.eE+-]+)", line)
            if m: results["obs"] = float(m.group(1))
        elif "Expected" in line and "%" in line:
            m = re.search(r"Expected\s+([\d\.]+).*?r\s*<\s*([0-9\.eE+-]+)", line)
            if m:
                pct = float(m.group(1))
                val = float(m.group(2))
                if pct == 2.5: results["0.025"] = val
                elif pct == 16.0: results["0.16"] = val
                elif pct == 50.0: results["0.5"] = val
                elif pct == 84.0: results["0.84"] = val
                elif pct == 97.5: results["0.975"] = val
    return results

def extract_name_and_mass(filename: str):
    m = re.match(r"datacard_(.+)_(\d+)\.txt$", filename)
    if m:
        name = f"{m.group(1)}_{m.group(2)}"
        return name, int(m.group(2))
    stem = Path(filename).stem
    return stem, None

def run_combine(datacard: Path, lumi: str, seed: str, logger):
    status_file = Path("status") / f"{datacard.stem}.Queue"
    if status_file.exists():
        status_file.rename(status_file.with_suffix(".Running"))
    
    logger.info(get_prefix() + f"Running {datacard.name} (AsymptoticLimits)")
    
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    
    r_max = estimate_rmax(datacard, safety_factor=10.0)
    
    cmd = [
        "combine",
        "-M", "AsymptoticLimits",
        datacard.name,
        "-s", str(seed),
        "-t", "-1", # Run expected limits (uses Asimov dataset)
        "--rMin", "0",
        "--rMax", str(r_max),
        "--setParameters", f"lumiscale={lumi}",
    ]
    
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / f"{datacard.stem}.log"
    
    try:
        with open(log_file, "w") as lf:
            res = subprocess.run(cmd, cwd=str(datacard.parent), env=env,
                                 stdout=lf, stderr=subprocess.STDOUT, text=True)
    except FileNotFoundError:
        raise SystemExit("combine not found in PATH. Activate your CMSSW env (cmsenv).")
        
    with open(log_file, "r") as lf:
        output = lf.read()
        
    limits = parse_limits(output)
    
    if limits:
        logger.info(get_prefix() + f"Complete {datacard.name} : median expected r < {limits.get('0.5', 'N/A')}")
    else:
        logger.info(get_prefix() + f"Failed {datacard.name} : return code {res.returncode}")
        
    run_file = status_file.with_suffix(".Running")
    if run_file.exists():
        run_file.rename(status_file.with_suffix(".Complete"))
        
    return limits, output, res.returncode, datacard.name

def cleanup(status_dir):
    # Asymptotic limits generate different root files than HybridNew
    for root_file in Path(".").glob("higgsCombine*.AsymptoticLimits.*.root"):
        root_file.unlink()

    log_dir = Path("logs")
    if log_dir.exists():
        merged_log = log_dir / "merged.log"
        with open(merged_log, "w") as mf:
            for log_file in sorted(log_dir.glob("*.log")):
                if log_file.name == "merged.log": continue
                with open(log_file, "r") as lf:
                    mf.write(f"--- {log_file.name} ---\n")
                    mf.write(lf.read())
                    mf.write("\n\n")
        logger.info(f"Merged logs into {merged_log}")
        
        # Remove individual log files
        for log_file in log_dir.glob("*.log"):
            if log_file != merged_log:
                log_file.unlink()
                
    if status_dir.exists():
        shutil.rmtree(status_dir)

def main():
    parser = argparse.ArgumentParser(description="Run combine AsymptoticLimits on datacards.")
    parser.add_argument("--lumi", default=DEFAULT_LUMI, help="Luminosity scale (default: 1)")
    parser.add_argument("--seed", default=DEFAULT_SEED, help="Random seed (default: 1)")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_JSON, help="Output JSON file (default: limits.json)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads (default: 4)")
    args = parser.parse_args()

    status_dir = Path("status")
    status_dir.mkdir(exist_ok=True)

    here = Path(".").resolve()
    txt_files = sorted(p for p in here.glob("*.txt"))
    if not txt_files:
        logger.info("No .txt files found in current directory.")
        return

    all_files = {}
    for f in txt_files:
        key, mass = extract_name_and_mass(f.name)
        all_files.setdefault(key, {})
        if mass is not None:
            all_files[key]["mass"] = mass

    # Prepare tasks (one per datacard, instead of per quantile)
    tasks = []
    for f in txt_files:
        tasks.append((f, args.lumi, args.seed))
        status_file = status_dir / f"{f.stem}.Queue"
        status_file.touch()

    # Run in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_task = {executor.submit(run_combine, datacard, lumi, seed, logger): datacard for datacard, lumi, seed in tasks}
        
        for future in concurrent.futures.as_completed(future_to_task):
            datacard = future_to_task[future]
            try:
                limits_dict, log, code, fname = future.result()
                key, _ = extract_name_and_mass(fname)
                if limits_dict:
                    all_files[key].update(limits_dict)
            except Exception as exc:
                logger.info(get_prefix() + f"Task for {datacard.name} generated an exception: {exc}")

    with open(args.output, "w") as fp:
        json.dump(all_files, fp, indent=2, sort_keys=True)
    logger.info(f"Wrote {args.output}")
    
    cleanup(status_dir)

if __name__ == "__main__":
    main()