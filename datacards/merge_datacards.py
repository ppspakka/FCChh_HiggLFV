#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Merge LFV Higgs datacards and submit limits.")
    parser.add_argument("--submit", action="store_true", help="Submit sbatch jobs after merging.")
    parser.add_argument("--njet-max", type=int, default=1, dest="njet_max", 
                        help="Maximum number of jets to merge (e.g. 1 merges 0j and 1j, 2 merges 0j, 1j, and 2j).")
    args = parser.parse_args()

    channels = ['etaumu', 'mutaue']
    njet_max = args.njet_max

    # Create an isolated environment for Combine scripts
    combine_env = os.environ.copy()
    combine_env["PYTHONNOUSERSITE"] = "1"

    for ch in channels:
        for cut in ["highmass", "lowmass"]:
            
            # Dynamically generate directory names from 0j to NJET_MAX
            jet_dirs = [f"datacards_{ch}_{cut}_{j}j" for j in range(njet_max + 1)]
            dir_comb = f"datacards_{ch}_{cut}_combined"

            # Check if all required input directories exist
            if not all(os.path.exists(d) for d in jet_dirs):
                print(f"Skipping {ch} {cut}: Missing one or more input directories (expected 0j to {njet_max}j).")
                continue

            os.makedirs(dir_comb, exist_ok=True)

            # Use 0j as the baseline to find all datacard filenames
            dir_0j = jet_dirs[0]
            cards_0j = glob.glob(os.path.join(dir_0j, "datacard_*.txt"))

            for card_0j in cards_0j:
                filename = os.path.basename(card_0j)
                card_comb = os.path.join(dir_comb, filename)

                # Gather paths for all jet multiplicities for this specific file
                card_paths = [card_0j]
                missing_card = False
                
                for j in range(1, njet_max + 1):
                    card_nj = os.path.join(jet_dirs[j], filename)
                    if not os.path.exists(card_nj):
                        print(f"Warning: {card_nj} missing. Skipping {filename}.")
                        missing_card = True
                        break
                    card_paths.append(card_nj)

                if missing_card:
                    continue

                # Execute combineCards.py dynamically mapping bin{j}j=path
                cmd = ["combineCards.py"]
                for j, cpath in enumerate(card_paths):
                    cmd.append(f"bin{j}j={cpath}")

                try:
                    # Capture the output as a string instead of writing directly to file
                    merged_output = subprocess.check_output(cmd, env=combine_env, universal_newlines=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running combineCards.py for {filename}")
                    continue

                # Fix the relative paths so they point to the parent directory ('../')
                fixed_lines = []
                for line in merged_output.splitlines():
                    if line.startswith("shapes "):
                        for d in jet_dirs:
                            line = line.replace(f"{d}/", f"../{d}/")
                    fixed_lines.append(line)

                # Write the modified output to the combined datacard
                with open(card_comb, "w") as f_out:
                    f_out.write("\n".join(fixed_lines) + "\n")

            print(f"Merged datacards created and paths adjusted in {dir_comb}/")

            # Copy auxiliary scripts from the 0j directory
            for script in ["run_limits.py", "slurm_submit.slurm"]:
                src = os.path.join(dir_0j, script)
                if os.path.exists(src):
                    shutil.copy(src, dir_comb)

            # Execute submission within the combined directory
            if args.submit:
                print(f"Submitting job for {ch} {cut} combined...")
                cwd = os.getcwd()
                os.chdir(dir_comb)
                # Submit using the standard environment, or combine_env if your Slurm script also invokes combine directly
                subprocess.run(["sbatch", "slurm_submit.slurm"], check=True)
                os.chdir(cwd)

if __name__ == "__main__":
    main()