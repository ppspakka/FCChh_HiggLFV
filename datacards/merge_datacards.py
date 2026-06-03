#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Merge LFV Higgs datacards and submit limits.")
    parser.add_argument("--submit", action="store_true", help="Submit sbatch jobs after merging.")
    args = parser.parse_args()

    channels = ['etaumu', 'mutaue', 'mue']

    # Create an isolated environment for Combine scripts
    combine_env = os.environ.copy()
    combine_env["PYTHONNOUSERSITE"] = "1"

    for ch in channels:
        for cut in ["highmass", "lowmass"]: # e.g datacards_etaumu_highmass_0j
            dir_0j = f"datacards_{ch}_{cut}_0j"
            dir_1j = f"datacards_{ch}_{cut}_1j"
            dir_comb = f"datacards_{ch}_{cut}_combined"

            if not os.path.exists(dir_0j) or not os.path.exists(dir_1j):
                print(f"Skipping {ch} {cut}: Missing input directories.")
                continue

            os.makedirs(dir_comb, exist_ok=True)

            cards_0j = glob.glob(os.path.join(dir_0j, "datacard_*.txt"))

            for card_0j in cards_0j:
                filename = os.path.basename(card_0j)
                card_1j = os.path.join(dir_1j, filename)
                card_comb = os.path.join(dir_comb, filename)

                if not os.path.exists(card_1j):
                    print(f"Warning: {card_1j} missing. Skipping {filename}.")
                    continue

                # Execute combineCards.py with isolated environment
                cmd = ["combineCards.py", f"bin0j={card_0j}", f"bin1j={card_1j}"]

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
                        line = line.replace(f"{dir_0j}/", f"../{dir_0j}/")
                        line = line.replace(f"{dir_1j}/", f"../{dir_1j}/")
                    fixed_lines.append(line)

                # Write the modified output to the combined datacard
                with open(card_comb, "w") as f_out:
                    f_out.write("\n".join(fixed_lines) + "\n")

            print(f"Merged datacards created and paths adjusted in {dir_comb}/")

            # Copy auxiliary scripts
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