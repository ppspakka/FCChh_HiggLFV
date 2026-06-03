#!/usr/bin/env python3
"""
Parallel runner for FCC-hh analysis

- Builds a queue of jobs across 8 pipeline configurations (2 channels * 2 mass categories * 2 jet bins).
- Supports multiple directories per signal/background (jobs run per directory, producing _idx files).
- Runs jobs in parallel.
- Automatically merges outputs via `hadd` and sums cutflow log statistics.
- Organizes or cleans up intermediate files.
- Loops makecard.py and slurm submission over all 8 final directories.
"""

from pathlib import Path
import subprocess
import sys
import os
import shutil
import concurrent.futures
import argparse
from dataclasses import dataclass
from typing import List, Dict, Any, Tuple
import re

# -------------------------
# User-configurable section
# -------------------------
PARALLEL = 4  # Default parallel workers

# Define the channels and categories (Total: 8 combinations)
CHANNELS = ["mutaue", "etaumu"]
# CATEGORIES = ["lowmass_0j", 
#             #   "lowmass_1j", "highmass_0j", "highmass_1j"
#               ]
# CHANNEL_CATS = [f"{ch}_{cat}" for ch in CHANNELS for cat in CATEGORIES]

NON_ORTHO_CUT = ["lowmass", 
                #  "highmass"
                 ]
NJET_MAX = 1
CHANNEL_CATS = []
for ch in CHANNELS:
    for cut in NON_ORTHO_CUT:
        for nj in range(NJET_MAX + 1):
            CHANNEL_CATS.append(f"{ch}_{cut}_{nj}j")

# Auto-generate pipeline paths based on naming convention
PIPELINES = {cc: f"./pipeline_{cc}.json" for cc in CHANNEL_CATS}

# Toggle which types to include
RUN_SIGNALS = True
RUN_BACKGROUNDS = True

# Signal Config (Modeled identically to background to support multi-dir structure)
SIGNAL_TYPE = "ggH"
MASS_RANGES = [200, 300, 450, 600, 750, 900]

# --- Replace these with your actual paths ---
SIGNALS = {
    f"{SIGNAL_TYPE}_{channel}_{mass}": [
        # f"/path/to/ggH/Hmass{mass}/DIR_1/",
        # f"/work/project/physics/psriling/FCC/FCChh/TestEnv/signals/{channel}/M450.root", # Test dir
        
        f"/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/signal/{channel}/M{mass}/",
    ]
    for mass in MASS_RANGES for channel in CHANNELS
}

BACKGROUNDS = {
    "ttbar": [
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/10Mseed10/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/10Mseed20/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/10Mseed30/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/10Mseed40/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/10Mseed50/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/5Mseed10/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/5Mseed20/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/5Mseed30/',
        # '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/5Mseed40/',
        '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/TThvq_leptonic/5Mseed50/',
    ],

    # "tW": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/ST_tW_topAndAnti-top_5f_inclusiveDecays_powheg/lep_mode/top_5Mseed10/',
    # ],

    # "tbarW": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/ST_tW_topAndAnti-top_5f_inclusiveDecays_powheg/lep_mode/antitop_5Mseed10/',
    # ],
    
    # "DY0j": [
    #     # Test dir
    #     # "/work/project/physics/psriling/FCC/FCChh/TestEnv/backgrounds/DY0j/dir1",
    #     # "/work/project/physics/psriling/FCC/FCChh/TestEnv/backgrounds/DY0j/dir2",
        
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/DY0Jets_tata/5Mseed100/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/DY0Jets_tata/5Mseed150/',
    # ],
    
    # "DY1j": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/DY1Jets_tata/5Mseed100/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/DY1Jets_tata/5Mseed150/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/DY1Jets_tata/5Mseed200/',
    # ],

    # "WW": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/WW_llvlvl/1Mseed100/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/WW_llvlvl/1Mseed110/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/WW_llvlvl/1Mseed120/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/WW_llvlvl/1Mseed130/',
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/WW_llvlvl/1Mseed140/',
    #     '/work/project/physics/psriling/FCC/Pythia8Delphes/Filter_ROOT/WW_llvlvl/5Mseed500/',
    #     '/work/project/physics/psriling/FCC/Pythia8Delphes/Filter_ROOT/WW_llvlvl/5Mseed600/',
    # ],
    
    # "SM_ggH_tautau": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/GluGluHToTauTau_powheg/1M_seed1/',
    # ],

    # "SM_ggH_WW": [
    #     '/work/project/physics/psriling/FCC/Pythia8Delphes/Filter_ROOT/GluGluHToWW_powheg/1M_seed1/',
    # ],

    # "SM_VBFH_tautau": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/VBF_H_powheg/tata/1M_seed1/',
    # ],

    # "SM_VBFH_WW": [
    #     '/work/project/cms/psriling/FCC_new/Pythia8Delphes/Filter_ROOT/VBF_H_powheg/WW/1M_seed2/',
    # ],
}

# --------------------------------------------

# Toggle per-category execution if needed (Defaults to globally True based on above)
CONFIG = {}
for cc in CHANNEL_CATS:
    CONFIG[f"signal_{cc}"] = RUN_SIGNALS
    CONFIG[f"background_{cc}"] = RUN_BACKGROUNDS


@dataclass
class Job:
    input_path: str
    idx: int
    output_dir: Path
    pipeline: str
    sample_type: str  # "signal" or "background"
    channel_cat: str  # e.g. "mutaue_lowmass_0j"
    name: str         # process name e.g. "ggH_200" or "DY0JET"

    def out_root(self) -> Path:
        return self.output_dir / f"{self.sample_type}_{self.name}_{self.idx}.root"

    def log_path(self) -> Path:
        return self.output_dir / f"{self.sample_type}_{self.name}_{self.idx}.log"

    def merged_root(self) -> Path:
        return self.output_dir / f"{self.sample_type}_{self.name}.root"

    def merged_log(self) -> Path:
        return self.output_dir / f"{self.sample_type}_{self.name}.log"

    def command(self) -> List[str]:
        # cpp_arg = f'analyze_pipeline.cpp("{self.input_path}","{self.out_root()}","{self.pipeline}")'
        cpp_arg = f'analyze_pipeline_optimize.cpp("{self.input_path}","{self.out_root()}","{self.pipeline}")'
        return ["root", "-l", "-b", "-q", cpp_arg]

# Setup pipeline by copying from current dir to output dir for reproducibility
# 1) copy the base pipeline files (ch_cut_0j.json)
# 2) Build the other jet pipeline by copying and replacing "0j" with e.g. "1j" in the filename and content
# 3) edit the content from line     "n_jet": 0, to    "n_jet": N for the respective pipeline
def setup_pipelines(base_out_dir: Path):
    # Setup base pipelines (0j)
    for ch in CHANNELS:
        for cut in NON_ORTHO_CUT:
            base_name = f"pipeline_{ch}_{cut}_0j.json"
            src_path = Path(base_name)
            if not src_path.exists():
                print(f"ERROR: Base pipeline file '{base_name}' not found.")
                continue
            dest_path = base_out_dir / src_path.name
            shutil.copy(src_path, dest_path)
            print(f"Copied base pipeline to {dest_path}")
    
    # Setup other jet pipelines by copying and modifying the base
    for ch in CHANNELS:
        for cut in NON_ORTHO_CUT:
            for nj in range(1, NJET_MAX + 1):
                base_name = f"pipeline_{ch}_{cut}_0j.json"
                new_name = f"pipeline_{ch}_{cut}_{nj}j.json"
                src_path = base_out_dir / base_name
                dest_path = base_out_dir / new_name
                if not src_path.exists():
                    print(f"ERROR: Base pipeline file '{src_path}' not found for creating '{new_name}'.")
                    continue
                shutil.copy(src_path, dest_path)
                
                # Modify the content to replace "0j" with "{nj}j" and update n_jet value
                with open(dest_path, "r") as f:
                    content = f.read()
                content = content.replace(f'"{base_name}"', f'"{new_name}"')
                content = re.sub(r'"n_jet":\s*0', f'"n_jet": {nj}', content)
                
                with open(dest_path, "w") as f:
                    f.write(content)
                
                print(f"Created pipeline {dest_path} with n_jet={nj}")



def update_job_status(status_dir: Path, job: Job, status: str):
    """Updates the status file for a job by replacing the old status file with a new one."""
    base_prefix = f"{job.channel_cat}_{job.name}_{job.idx}__."
    for f in status_dir.glob(f"{base_prefix}*"):
        f.unlink(missing_ok=True)
    (status_dir / f"{base_prefix}{status}").touch()


def build_jobs(base_out_dir: Path) -> List[Job]:
    jobs: List[Job] = []
    
    for cc in CHANNEL_CATS:
        out_dir = base_out_dir / f"hist_{cc}"
        
        # Signals
        # if CONFIG.get(f"signal_{cc}", False):
        #     for sig_name, paths in SIGNALS.items():
        #         for idx, path in enumerate(paths):
        #             jobs.append(Job(
        #                 input_path=path, idx=idx, output_dir=out_dir, 
        #                 pipeline=PIPELINES[cc], sample_type="signal", 
        #                 channel_cat=cc, name=sig_name
        #             ))
        
        # Only run signal mutaue for channels containing "mutaue"
        if CONFIG.get(f"signal_{cc}", False) and "mutaue" in cc:
            for sig_name, paths in SIGNALS.items():
                # if mutaue in sig_name
                if sig_name.startswith(f"{SIGNAL_TYPE}_mutaue"):
                    for idx, path in enumerate(paths):
                        jobs.append(Job(
                            input_path=path, idx=idx, output_dir=out_dir, 
                            # pipeline=PIPELINES[cc], sample_type="signal", 
                            # pipeline in path base_out_dir/PIPELINES[cc] via path.join
                            pipeline=os.path.join(base_out_dir, PIPELINES[cc]), sample_type="signal",
                            channel_cat=cc, name=sig_name
                        ))
        # Etaumu
        if CONFIG.get(f"signal_{cc}", False) and "etaumu" in cc:
            for sig_name, paths in SIGNALS.items():
                if sig_name.startswith(f"{SIGNAL_TYPE}_etaumu"):
                    for idx, path in enumerate(paths):
                        jobs.append(Job(
                            input_path=path, idx=idx, output_dir=out_dir, 
                            pipeline=os.path.join(base_out_dir, PIPELINES[cc]), sample_type="signal", 
                            channel_cat=cc, name=sig_name
                        ))

        # Backgrounds
        if CONFIG.get(f"background_{cc}", False):
            for bg_name, paths in BACKGROUNDS.items():
                for idx, path in enumerate(paths):
                    jobs.append(Job(
                        input_path=path, idx=idx, output_dir=out_dir, 
                        pipeline=os.path.join(base_out_dir, PIPELINES[cc]), sample_type="background", 
                        channel_cat=cc, name=bg_name
                    ))
    return jobs


def run_job(job: Job, status_dir: Path) -> Tuple[Job, int, str]:
    """Runs a job, captures its output, and manages its status."""
    update_job_status(status_dir, job, "running")
    job.output_dir.mkdir(parents=True, exist_ok=True)
    
    proc = subprocess.run(job.command(), capture_output=True, text=True)
    
    with open(job.log_path(), "w", encoding="utf-8", errors="replace") as logf:
        logf.write(proc.stdout)
        if proc.stderr:
            logf.write("\n--- STDERR ---\n")
            logf.write(proc.stderr)

    update_job_status(status_dir, job, "done")
    return job, proc.returncode, proc.stdout


def merge_and_cleanup(jobs: List[Job], base_out_dir: Path, clean_intermediates: bool):
    """Merges root/log files from multiple directories and handles intermediates."""
    print("\n" + "="*50)
    print("Starting Merging & Cleanup Process")
    print("="*50)
    
    # Group jobs by logical output (channel_cat, sample_type, name)
    grouped_jobs: Dict[str, List[Job]] = {}
    for j in jobs:
        key = f"{j.channel_cat}|{j.sample_type}|{j.name}"
        if key not in grouped_jobs:
            grouped_jobs[key] = []
        grouped_jobs[key].append(j)

    for key, group in grouped_jobs.items():
        if not group:
            continue
            
        first_job = group[0]
        channel_cat = first_job.channel_cat
        merged_root = first_job.merged_root()
        merged_log = first_job.merged_log()
        
        # 1. Merge ROOT files using hadd
        root_files = [str(j.out_root()) for j in group if j.out_root().exists()]
        if not root_files:
            print(f"Warning: No ROOT files found to merge for {key}")
            continue
            
        if len(root_files) == 1:
            # Only one file, just rename/copy
            shutil.copy(root_files[0], str(merged_root))
            print(f"Copied single root file for {first_job.name} ({channel_cat})")
        else:
            # Multiple files, hadd required
            cmd = ["hadd", "-f", str(merged_root)] + root_files
            print(f"Merging {len(root_files)} ROOT files for {first_job.name} ({channel_cat})")
            subprocess.run(cmd, capture_output=True, check=True)

        # 2. Merge Log files & Sum Cutflows
        total_events = 0
        cut_counts = {}
        cut_order = []
        
        # Regex mappings based on expected output log format
        re_total = re.compile(r"Total events:\s+(\d+)")
        re_cut = re.compile(r"After\s+(.*?)\s+(\d+)\s+\(Efficiency")
        
        for j in group:
            if not j.log_path().exists():
                continue
            with open(j.log_path(), "r") as f:
                for line in f:
                    m_total = re_total.search(line)
                    if m_total:
                        total_events += int(m_total.group(1))
                    
                    m_cut = re_cut.search(line)
                    if m_cut:
                        c_name = m_cut.group(1).strip()
                        c_val = int(m_cut.group(2))
                        if c_name not in cut_order:
                            cut_order.append(c_name)
                        cut_counts[c_name] = cut_counts.get(c_name, 0) + c_val

        # Write merged log
        with open(merged_log, "w") as f:
            f.write("==== Merged Analysis Parameters ====\n")
            f.write(f"(Merged across {len(group)} input directories)\n")
            f.write("======================================\n\n")
            f.write("==== Pipeline summary (Merged) ====\n")
            f.write(f"Total events:          {total_events}\n")
            
            for cut in cut_order:
                val = cut_counts[cut]
                eff = (val / total_events * 100.0) if total_events > 0 else 0.0
                f.write(f"After {cut:<25} {val} (Efficiency: {eff:.6f}%)\n")

        # 3. Handle Intermediates (Delete or Stash)
        inter_dir = base_out_dir / "intermediates" / channel_cat
        if not clean_intermediates:
            inter_dir.mkdir(parents=True, exist_ok=True)
            
        for j in group:
            if j.out_root().exists():
                if clean_intermediates:
                    j.out_root().unlink()
                else:
                    shutil.move(str(j.out_root()), str(inter_dir / j.out_root().name))
            
            if j.log_path().exists():
                if clean_intermediates:
                    j.log_path().unlink()
                else:
                    shutil.move(str(j.log_path()), str(inter_dir / j.log_path().name))

    print("Merging & Cleanup Complete.")


def run_makecard_commands(args, dry_run: bool = False):
    print("\n" + "="*30)
    print("Starting makecard generation")
    print("="*30)

    lumi_pb = 30_000_000  # FCC-hh expected lumi scale? Adjust as needed
    failures = 0
    base_out = Path(args.output_dir)

    for cc in CHANNEL_CATS:
        in_dir = base_out / f"hist_{cc}"
        out_dir = base_out / f"datacards_{cc}"
        
        # cut_type = "lowmass" if "lowmass" in cc else "highmass"
        cut_type = "lowmass"
        
        # Skip if input dir wasn't created or is empty
        if not in_dir.exists() or not any(in_dir.iterdir()):
            continue
            
        out_dir.mkdir(exist_ok=True)
        
        cmd = [
            "python3", "makecard.py",
            "--in-dir", str(in_dir),
            "--lumi-pb", str(lumi_pb),
            "--out-root", str(out_dir / "merged.root"),
            "--out-card", str(out_dir),
            "--bin-type", cut_type
        ]

        print(f"\nRunning makecard for '{cc}':")
        if dry_run:
            print("DRY RUN CMD:", " ".join(cmd))
            continue

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Makecard for '{cc}' completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: makecard for '{cc}' failed (rc {e.returncode}).")
            print("STDERR:", e.stderr)
            failures += 1
        except FileNotFoundError:
            print("ERROR: 'makecard.py' not found.")
            failures += 1
            break
    
    if failures > 0:
        print(f"\n{failures} makecard job(s) failed.")
        sys.exit(3)


def run_sbatch_commands(args):
    script_files = ["datacards/run_limits.py", "datacards/slurm_submit.slurm"]
    base_out = Path(args.output_dir)
    
    for cc in CHANNEL_CATS:
        in_dir = base_out / f"hist_{cc}"
        out_dir = base_out / f"datacards_{cc}"
        
        # Only copy if datacard dir was actually created
        if not out_dir.exists():
            continue
            
        for script in script_files:
            if not Path(script).exists():
                continue
            dest = out_dir / Path(script).name
            try:
                shutil.copy(script, dest)
                print(f"Copied {script} to {dest}")
            except Exception as e:
                print(f"Failed to copy {script} to {dest}: {e}")
    


def main():
    parser = argparse.ArgumentParser(description="Parallel FCC-hh Runner")
    parser.add_argument("--output-dir", "-o", type=str, default="./fcc_hh_cuts", help="parent output directory")
    parser.add_argument("--parallel", "-p", type=int, default=PARALLEL, help="number of parallel jobs")
    parser.add_argument("--list", action="store_true", help="only list jobs (don't execute)")
    parser.add_argument("--dry-run", action="store_true", help="show commands without running")
    parser.add_argument("--clean-intermediates", action="store_true", help="delete _idx root/log files instead of stashing them")
    parser.add_argument("--skip-makecard", action="store_true", help="skip makecard step")
    parser.add_argument("--skip-sbatch", action="store_true", help="skip sbatch submission step")
    parser.add_argument("--skip-run", action="store_true", help="skip the job execution step (useful for just merging/makecards)")
    args = parser.parse_args()

    out_dir_path = Path(args.output_dir)
    status_dir = out_dir_path / "status"
    
    jobs = build_jobs(out_dir_path)

    if not jobs:
        print("No jobs built. Check configuration.")
        return

    print(f"Built {len(jobs)} jobs. Outputting to {args.output_dir}")

    if args.list:
        for j in jobs:
            print(f"CMD: {' '.join(j.command())} -> {j.out_root()}")
        return

    if args.dry_run:
        for j in jobs:
            print(f"DRY RUN: {' '.join(j.command())}")
        if not args.skip_makecard:
            run_makecard_commands(args, dry_run=True)
        return

    failures = []
    if not args.skip_run:
        if status_dir.exists():
            shutil.rmtree(status_dir)
        status_dir.mkdir(parents=True, exist_ok=True)
        
        # # Copy pipelines to output for reproducibility
        # for cc, pf in PIPELINES.items():
        #     if Path(pf).exists():
        #         shutil.copy(pf, out_dir_path / Path(pf).name)
        setup_pipelines(out_dir_path)
        
        for job in jobs:
            update_job_status(status_dir, job, "queue")

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as exe:
            future_to_job = {exe.submit(run_job, job, status_dir): job for job in jobs}
            
            for fut in concurrent.futures.as_completed(future_to_job):
                job = future_to_job[fut]
                try:
                    completed_job, rc, out_text = fut.result()
                    
                    n_queue = len(list(status_dir.glob("*.queue")))
                    n_run = len(list(status_dir.glob("*.running")))
                    n_done = len(list(status_dir.glob("*.done")))
                    
                    print(f"\n{'='*75}")
                    print(f"[Queue {n_queue}] [Running {n_run}] [Done {n_done}] --- Log for: {job.channel_cat} | {job.name}_{job.idx}")
                    print(f"{'='*75}")
                    print(out_text.strip())
                    
                    if rc != 0:
                        print(f"-> FAILED (rc={rc}): {job.input_path}")
                        failures.append((job, rc))
                    else:
                        print(f"-> DONE: {job.input_path}")

                except Exception as e:
                    print(f"EXCEPTION for {job.input_path}: {e}")
                    failures.append((job, -1))
        

    if failures:
        print(f"\n{len(failures)} jobs failed. Check logs. Aborting merge and downstream steps.")
        sys.exit(2)
        
    # Execute Merge & Cleanup logic
    merge_and_cleanup(jobs, out_dir_path, args.clean_intermediates)

    # Downstream processes
    if not args.skip_makecard:
        run_makecard_commands(args)
    
    run_sbatch_commands(args)
    
    # Run Final merge_datacards script
    merge_script = "datacards/merge_datacards.py"
    if Path(merge_script).exists():
        dest_merge = out_dir_path / Path(merge_script).name
        shutil.copy(merge_script, dest_merge)
        
        cmd = ["python3", "merge_datacards.py", "--njet-max", str(NJET_MAX)]
        if not args.skip_sbatch:
            cmd.append("--submit")
            
        print("\nRunning merge_datacards.py with command:", " ".join(cmd))
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=out_dir_path)
            print("merge_datacards.py completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: merge_datacards.py failed (rc={e.returncode}).")
            print(e.stderr)
            sys.exit(4)

if __name__ == "__main__":
    main()