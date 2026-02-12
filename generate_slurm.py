#!/usr/bin/env python3
"""Interactive SLURM script generator for anchored docking with tab completion."""

import glob
import os
import readline
import sys


def path_completer(text, state):
    """Tab-complete file/directory paths."""
    if "~" in text:
        text = os.path.expanduser(text)
    if text == "":
        matches = glob.glob("*")
    else:
        matches = glob.glob(text + "*")
    matches = [m + "/" if os.path.isdir(m) else m for m in matches]
    try:
        return matches[state]
    except IndexError:
        return None


def choice_completer(choices):
    """Return a completer function for a fixed set of choices."""
    def completer(text, state):
        matches = [c for c in choices if c.startswith(text)]
        try:
            return matches[state]
        except IndexError:
            return None
    return completer


def set_completer(func):
    readline.set_completer(func)
    readline.set_completer_delims(" \t\n")
    readline.parse_and_bind("tab: complete")


def prompt_path(label, default=None):
    set_completer(path_completer)
    hint = f" [{default}]" if default else ""
    val = input(f"  {label}{hint}: ").strip()
    return val if val else default


def prompt_choice(label, choices, default=None):
    set_completer(choice_completer(choices))
    hint = f" [{default}]" if default else ""
    while True:
        val = input(f"  {label} ({'/'.join(choices)}){hint}: ").strip()
        val = val if val else default
        if val in choices:
            return val
        print(f"    Invalid choice. Options: {', '.join(choices)}")


def prompt_value(label, default=None, cast=None):
    readline.set_completer(None)
    hint = f" [{default}]" if default else ""
    val = input(f"  {label}{hint}: ").strip()
    val = val if val else default
    if val is not None and cast is not None:
        try:
            val = cast(val)
        except (ValueError, TypeError):
            print(f"    Invalid value, using default: {default}")
            val = default
    return val


def prompt_bool(label, default=False):
    set_completer(choice_completer(["yes", "no"]))
    hint = "yes" if default else "no"
    val = input(f"  {label} (yes/no) [{hint}]: ").strip().lower()
    if val in ("yes", "y"):
        return True
    if val in ("no", "n"):
        return False
    return default


def main():
    print("\n=== Anchored Docking SLURM Script Generator ===\n")
    print("Press Tab to autocomplete paths and choices.\n")

    # --- SLURM settings ---
    print("-- SLURM settings --")
    partition = prompt_value("Partition", default="ga100")
    gpu_count = prompt_value("Number of GPUs", default="1", cast=int)
    job_name = prompt_value("Job name", default="anchored_dock")
    time_limit = prompt_value("Time limit (e.g. 2:00:00)", default=None)
    modules = prompt_value(
        "Modules to load (comma-separated)",
        default="pixi/0.56.0,AutoDock-GPU/1.5.3-CUDA",
    )

    # --- Required docking parameters ---
    print("\n-- Docking parameters (required) --")
    template_file = prompt_path("Template ligand SDF (--template_file)")
    print("\n  This can be a single molecule or a multi-molecule SDF file.")
    print("  If it contains multiple molecules, each one will be docked separately.")
    input_file = prompt_path("Input ligand SDF (--input_file)")
    receptor_file = prompt_path("Receptor PDB (--receptor_file)")
    docking_dir = prompt_path("Output directory (--docking_dir)", default="output")

    # --- Paths ---
    print("\n-- Tool paths --")
    adg_path = prompt_path("AutoDock-GPU bin dir (--adg_path)", default="$PWD/bin")
    adfr_path = prompt_path("ADFRsuite dir (--adfr_path)", default="/nemo/stp/chemicalbiology/home/shared/software/ADFRsuite")
    scripts_path = prompt_path("Scripts dir (--scripts_path)", default="scripts")

    # --- Optional docking parameters ---
    print("\n-- Optional docking parameters (press Enter to keep defaults) --")

    print("\n  Anchors are atoms shared between your template and input molecule")
    print("  that act like pins to hold the molecule in the right orientation.")
    anchor_num = prompt_value("How many anchor pins to use", default="3", cast=int)

    print("\n  'random' picks anchor atoms randomly from the shared structure.")
    print("  'area' picks 3 anchors that are as spread out as possible (like a tripod).")
    anchor_mode = prompt_choice("Anchor selection strategy", ["random", "area"], default="random")

    print("\n  The search box is the region around the binding site where the")
    print("  molecule is allowed to move. Larger = more freedom, slower search.")
    box_size = prompt_value("Search box size in Angstroms", default="20.0", cast=float)

    print("\n  Normally the shared structure between template and input is found")
    print("  automatically. Set this only if you want to manually define which")
    print("  part of the molecule to match (advanced, uses SMARTS notation).")
    query_smarts = prompt_value("Manual substructure pattern (SMARTS)", default=None)

    if query_smarts:
        print("\n  Since you set a manual SMARTS, you can pick exactly which atoms")
        print("  in that pattern to use as anchors (0-based atom positions).")
        anchor_indices = prompt_value("Specific anchor atom positions (space-separated)", default=None)
    else:
        anchor_indices = None

    if anchor_mode == "area":
        print("\n  When using the 'area' strategy, the script tries many random sets")
        print("  of 3 anchors and keeps the most spread-out set. More samples = better")
        print("  result but slower.")
        samples = prompt_value("Number of random attempts for tripod selection", default="100", cast=int)
    else:
        samples = 100

    print("\n  Set a seed to get the same anchor selection every time you run.")
    print("  Leave blank for a different random pick each run.")
    random_seed = prompt_value("Random seed for reproducibility", default=None)

    force = prompt_bool("Overwrite output directory if it exists", default=True)
    clean = prompt_bool("Delete intermediate files after docking finishes", default=False)
    verbose = prompt_bool("Show extra debug info during the run", default=False)

    # --- Build the command ---
    cmd_parts = [
        "pixi run python anchored_docking.py",
        f"    --template_file {template_file}",
        f"    --input_file {input_file}",
        f"    --receptor_file {receptor_file}",
        f"    --docking_dir {docking_dir}",
        f"    --adg_path {adg_path}",
        f"    --adfr_path {adfr_path}",
        f"    --scripts_path {scripts_path}",
        f"    --anchor_num {anchor_num}",
        f"    --anchor_mode {anchor_mode}",
        f"    --size {box_size}",
        f"    --samples {samples}",
    ]
    if query_smarts:
        cmd_parts.append(f"    --query_smarts '{query_smarts}'")
    if anchor_indices:
        cmd_parts.append(f"    --anchor_indices {anchor_indices}")
    if random_seed is not None:
        cmd_parts.append(f"    --random_seed {random_seed}")
    if force:
        cmd_parts.append("    --force")
    if clean:
        cmd_parts.append("    --clean")
    if verbose:
        cmd_parts.append("    --verbose")

    dock_cmd = " \\\n".join(cmd_parts)

    # --- Build the SLURM script ---
    lines = [
        "#!/bin/bash",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --gres=gpu:{gpu_count}",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --output=dock_%j.log",
    ]
    if time_limit:
        lines.append(f"#SBATCH --time={time_limit}")
    lines.append("")
    for mod in modules.split(","):
        mod = mod.strip()
        if mod:
            lines.append(f"module load {mod}")
    lines.append("")
    lines.append(dock_cmd)
    lines.append("")

    script_content = "\n".join(lines)

    # --- Preview and save ---
    print("\n--- Generated SLURM script ---")
    print(script_content)

    set_completer(path_completer)
    output_file = prompt_path("Save as", default="submit_dock.sh")
    if output_file:
        with open(output_file, "w") as f:
            f.write(script_content)
        os.chmod(output_file, 0o755)
        print(f"\nSaved to {output_file}")
        print(f"Submit with: sbatch {output_file}")
    else:
        print("\nScript not saved.")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print("\n\nCancelled.")
        sys.exit(0)
