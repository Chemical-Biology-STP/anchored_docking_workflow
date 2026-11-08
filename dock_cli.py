#!/usr/bin/env python3
"""Interactive command-line interface for the anchored docking workflow.

Start it with:  pixi run dock

It walks through the inputs one at a time, checks each of them, then either
submits the docking job to SLURM, runs it in the current terminal, or just
saves the job script. It can also watch a submitted job until it finishes.
"""

import argparse
import getpass
import glob
import math
import os
import re
import readline
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_DIR)

# The lab's launcher (ddsuite) sets DDSUITE_RUN_COMMAND, e.g. "ddsuite dock". It is used in messages that tell the
# user how to start this tool; on its own the tool is started with `pixi run dock`.
UNDER_DDSUITE = bool(os.environ.get("DDSUITE_RUN_COMMAND"))
RUN_COMMAND = os.environ.get("DDSUITE_RUN_COMMAND") or "pixi run dock"

DEFAULT_MODULES = "GCC/13.2.0,AutoDock-GPU/1.5.3-CUDA"
DEFAULT_PARTITION = "ga100"
# CUDA compute capability of the GPUs behind each known partition; a custom AutoDock-GPU build only runs on
# the capabilities it was built for (its `targets` file). A partition not listed here is not checked.
PARTITION_COMPUTE = {"ga100": "80", "gh100": "90"}
# Partitions that can only be used together with a reservation.
PARTITION_RESERVATION = {"gh100": "h100"}
SUBMIT_CHECK_SECONDS = 30   # how long to wait for Slurm when checking access
# Deliberately generous: the Bace1 example takes about 2 minutes per ligand.
MINUTES_PER_LIGAND = 10
MIN_TIME_MINUTES = 30
MIN_CORE_ATOMS = 6  # smaller than one ring: unlikely to be the same chemical series
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "TIP3", "SOL"}
TOTAL_STEPS = 6

# Limits compiled into the AutoDock-GPU 1.5.3 build of the cluster module (common/defines.h).
MAX_ATOM_TYPES = 32        # MAX_NUM_OF_ATYPES; scripts/build_autodock_gpu.sh makes a build with a higher value
CUSTOM_BUILD_GLOB = "adgpu_atypes*"   # folders in <repo>/bin written by that script
# Where the group's shared AutoDock-GPU installations live; a build with a higher atom-type limit installed
# there (module AutoDock-GPU/1.5.3-CUDA-atypes64) is found as <root>/<version>/bin.
SHARED_BUILD_ROOT = "/nemo/stp/chemicalbiology/home/shared/easybuild/software/AutoDock-GPU"
MAX_GRID_POINTS = 256      # MAX_NUM_GRIDPOINTS
# Largest box (Angstrom) for which grid_points() stays within MAX_GRID_POINTS.
MAX_BOX_SIZE = 95.0
# Grid generation and disk use grow with the cube of the box size.
LARGE_BOX_SIZE = 40.0
# Grid maps and other files for one ligand with the default 20 A box (measured on the example).
MB_PER_LIGAND_20A = 30.0
MB_PER_ANCHOR_20A = 1.3    # each anchor adds a grid map: 24 MB with 3 anchors, 41 MB with 17 on the example
MB_KEPT_PER_LIGAND = 1.0   # files kept per ligand when intermediates are deleted

POLL_SECONDS = 10
MAX_WATCH_ERRORS = 6       # consecutive failed squeue calls before giving up watching
ACTIVE_STATES = {"PENDING", "RUNNING", "REQUEUED", "RESIZING", "SUSPENDED", "COMPLETING", "CONFIGURING"}

# --------------------------------------------------------------------------
# Terminal output
# --------------------------------------------------------------------------

USE_COLOR = sys.stdout.isatty() and "NO_COLOR" not in os.environ
WARNINGS = []


def _c(code, text):
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else text


def bold(text):
    return _c("1", text)


def dim(text):
    return _c("2", text)


def say(msg=""):
    print(msg)


def ok(msg):
    print(f"  {_c('32', '[ok]')} {msg}")


def warn(msg, record=True):
    """Print a warning. Recorded warnings are repeated in the final review."""
    if record:
        WARNINGS.append(msg)
    print(f"  {_c('33', '[!]')}  {msg}")


def fail(msg):
    print(f"  {_c('31', '[x]')}  {msg}")


def note(msg):
    print(f"       {dim(msg)}")


def heading(step, title):
    text = f"Step {step} of {TOTAL_STEPS}: {title}"
    print(f"\n{bold(text)}\n{'-' * len(text)}")


def table(rows, headers=None, indent="  "):
    rows = [[str(c) for c in r] for r in rows]
    if headers:
        rows = [list(headers)] + rows
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for n, r in enumerate(rows):
        print(indent + "  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip())
        if headers and n == 0:
            print(indent + "  ".join("-" * w for w in widths))


# --------------------------------------------------------------------------
# Prompts
# --------------------------------------------------------------------------


def _set_completer(func):
    readline.set_completer(func)
    readline.set_completer_delims(" \t\n")
    readline.parse_and_bind("tab: complete")


def path_completer(text, state):
    text = os.path.expanduser(text)
    matches = glob.glob(text + "*")
    matches = [m + "/" if os.path.isdir(m) else m for m in matches]
    try:
        return sorted(matches)[state]
    except IndexError:
        return None


def choice_completer(choices):
    def completer(text, state):
        matches = [c for c in choices if c.startswith(text)]
        try:
            return matches[state]
        except IndexError:
            return None

    return completer


def ask(label, default=None, validate=None, completer=None, required=True):
    """Prompt until the answer is valid. validate() may raise ValueError."""
    _set_completer(completer)
    suffix = f" [{default}]" if default not in (None, "") else ""
    while True:
        raw = input(f"  {label}{suffix}: ").strip()
        value = raw if raw else default
        if value in (None, ""):
            if required:
                print("       Please enter a value.")
                continue
            return None
        if validate is not None:
            try:
                value = validate(value)
            except ValueError as exc:
                fail(str(exc))
                continue
        return value


def ask_yes_no(label, default=False):
    def validate(text):
        t = str(text).lower()
        if t in ("y", "yes", "true"):
            return True
        if t in ("n", "no", "false"):
            return False
        raise ValueError("Please answer yes or no.")

    return ask(f"{label} (yes/no)", "yes" if default else "no", validate, choice_completer(["yes", "no"]))


def ask_choice(label, choices, default):
    def validate(text):
        if text not in choices:
            raise ValueError(f"Please choose one of: {', '.join(choices)}")
        return text

    return ask(f"{label} ({'/'.join(choices)})", default, validate, choice_completer(choices))


def ask_number(label, default, cast=float, minimum=None, required=True, maximum=None):
    def validate(text):
        try:
            value = cast(text)
        except (TypeError, ValueError):
            raise ValueError("Please enter a number.")
        if minimum is not None and value < minimum:
            raise ValueError(f"Please enter a number of at least {minimum}.")
        if maximum is not None and value > maximum:
            raise ValueError(f"Please enter a number of at most {maximum}.")
        return value

    return ask(label, None if default is None else str(default), validate, required=required)


def file_validator(extensions=()):
    def validate(text):
        path = os.path.abspath(os.path.expanduser(os.path.expandvars(text.strip("'\""))))
        if not os.path.isfile(path):
            raise ValueError(f"File not found: {path}")
        if extensions and not path.lower().endswith(extensions):
            raise ValueError(f"Expected a {' or '.join(extensions)} file.")
        if os.path.getsize(path) == 0:
            raise ValueError("That file is empty.")
        return path

    return validate


def time_validator(text):
    if not re.fullmatch(r"(\d+-)?\d+(:\d{2}){0,2}", text):
        raise ValueError("Use hours:minutes:seconds, for example 2:00:00 (or days-hours:minutes:seconds).")
    return text


def fmt_minutes(minutes):
    minutes = int(math.ceil(minutes))
    days, rest = divmod(minutes, 1440)
    hours, mins = divmod(rest, 60)
    return f"{days}-{hours:02d}:{mins:02d}:00" if days else f"{hours}:{mins:02d}:00"


# --------------------------------------------------------------------------
# Limits, resources and environment
# --------------------------------------------------------------------------


def grid_points(size):
    """Grid points per axis for a box of `size` A, as scripts/write-gpf.py sets it up (npts + 1)."""
    return 2 * int(size / 0.75) + 1


def find_custom_build():
    """The AutoDock-GPU build with the highest atom-type limit made by scripts/build_autodock_gpu.sh, or None.

    Looks in the group's shared software folder (a cluster module) and in <repo>/bin.

    Returns (folder, limit).
    """
    import anchored_docking as ad

    best = None
    folders = sorted(glob.glob(os.path.join(SHARED_BUILD_ROOT, "*", "bin"))) \
        + sorted(glob.glob(os.path.join(REPO_DIR, "bin", CUSTOM_BUILD_GLOB)))   # the shared one wins a tie
    for folder in folders:
        binary = next((os.path.join(folder, n) for n in ad.ADGPU_BINARY_NAMES
                       if os.path.isfile(os.path.join(folder, n)) and os.access(os.path.join(folder, n), os.X_OK)), None)
        if binary is None:
            continue
        limit = ad.atom_type_limit(binary)
        if limit > MAX_ATOM_TYPES and (best is None or limit > best[1]):
            best = (folder, limit)
    return best


def build_targets(folder):
    """Compute capabilities a custom build was made for (from its `targets` file), or None if unknown."""
    try:
        with open(os.path.join(folder, "targets")) as fh:
            return fh.read().split()
    except OSError:
        return None


@dataclass
class SubmitCheck:
    ok: object          # True, False, or None when it could not be found out
    reason: str = ""


# Replies of the scheduler that mean "try again later" and say nothing about access.
SCHEDULER_UNAVAILABLE = ("unable to contact", "socket timed out", "connection refused", "slurm controller",
                         "temporarily")


def explain_submit_error(text, partition, reservation=""):
    """Say in plain words why Slurm refused a job, from its error text."""
    low = text.lower()
    if "access denied to requested reservation" in low:
        return f"you are not allowed to use the reservation {reservation}"
    if "requested reservation is invalid" in low:
        return f"the reservation {reservation} does not exist"
    if "invalid partition" in low:
        return f"the partition {partition} does not exist"
    if "invalid account" in low or "account/partition" in low:
        return f"your account cannot use the partition {partition}"
    if "invalid qos" in low or "qos" in low and "not permitted" in low:
        return f"your account has no access to the quality of service that the partition {partition} needs"
    if "time limit" in low:
        return f"the time limit is longer than the partition {partition} allows"
    if "node configuration is not available" in low:
        needed = PARTITION_RESERVATION.get(partition)
        if needed and not reservation:
            return f"the partition {partition} can only be used with the reservation {needed}"
        return f"no node of the partition {partition} has the requested GPUs"
    lines = [re.sub(r"^(sbatch: )?(error|allocation failure): ?", "", l.strip(), flags=re.I)
             for l in text.splitlines() if l.strip()]
    return lines[0] if lines else "the scheduler refused the job"


def check_submission(partition, reservation="", gpus=1, time_limit="0:10:00"):
    """Ask Slurm whether this account could submit such a job (sbatch --test-only: nothing is queued).

    Access to a partition can depend on the account and on reservations, so Slurm's own rules are used instead of
    a list kept here.
    """
    sbatch = shutil.which("sbatch")  # the real program; the interactive shell wrapper is not visible from here
    if sbatch is None:
        return SubmitCheck(None, "sbatch was not found")
    command = [sbatch, "--test-only", f"--partition={partition}", f"--gres=gpu:{gpus}", f"--time={time_limit}",
               "--wrap=true"]
    if reservation:
        command.insert(2, f"--reservation={reservation}")
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=SUBMIT_CHECK_SECONDS,
                                stdin=subprocess.DEVNULL, check=False)
    except (subprocess.TimeoutExpired, OSError):
        return SubmitCheck(None, "the scheduler did not answer")
    if result.returncode == 0:
        return SubmitCheck(True)
    text = result.stdout + result.stderr
    if any(word in text.lower() for word in SCHEDULER_UNAVAILABLE):
        return SubmitCheck(None, "the scheduler could not be reached")
    return SubmitCheck(False, explain_submit_error(text, partition, reservation))


def usable_gpu_partitions(s, exclude=()):
    """The known GPU partitions this account can submit to, according to Slurm."""
    found = []
    for name in PARTITION_COMPUTE:
        if name in exclude:
            continue
        if check_submission(name, PARTITION_RESERVATION.get(name, ""), s.gpus, s.time_limit).ok is True:
            found.append(name)
    return found


def ask_reservation(s):
    return ask("Reservation (type none for no reservation)", PARTITION_RESERVATION.get(s.partition, "none"),
               lambda t: "" if t.lower() == "none" else t)


def ensure_partition_access(s):
    """Before the job script is written, make sure this account can use the chosen partition and reservation.

    Not every account may use every partition (the H100 nodes need a reservation that only some users are on).
    """
    if s.mode == "run" or shutil.which("sbatch") is None:
        return
    while True:
        result = check_submission(s.partition, s.reservation, s.gpus, s.time_limit)
        if result.ok is True:
            return
        if result.ok is None:
            note(f"Could not check whether you may use the partition {s.partition}: {result.reason}.")
            return
        with_res = f" with the reservation {s.reservation}" if s.reservation else ""
        fail(f"You cannot submit to the partition {s.partition}{with_res}: {result.reason}.")
        usual = PARTITION_RESERVATION.get(s.partition, "")
        if usual != s.reservation and check_submission(s.partition, usual, s.gpus, s.time_limit).ok is True:
            say(f"  The partition {s.partition} works for you with {'the reservation ' + usual if usual else 'no reservation'}.")
        others = usable_gpu_partitions(s, exclude=(s.partition,))
        if others:
            say(f"  You can use: {', '.join(others)}.")
        if ask_yes_no("Use it anyway (the job will be refused when it is submitted)?", False):
            warn(f"Slurm will refuse jobs for the partition {s.partition}{with_res}: {result.reason}.")
            return
        s.partition = ask("Partition", others[0] if others else DEFAULT_PARTITION)
        s.reservation = ask_reservation(s)


def build_problem(s):
    """Why the chosen custom AutoDock-GPU would not run on the chosen partition, or None."""
    if not s.adg_path:
        return None
    wanted = PARTITION_COMPUTE.get(s.partition)
    have = build_targets(s.adg_path)
    if wanted is None or have is None:
        return None
    if wanted not in have:
        return (f"The custom AutoDock-GPU in {s.adg_path} was built for compute capability {' '.join(have)}, but "
                f"partition {s.partition} has GPUs of capability {wanted}. It would fail there.")
    return None


def atom_type_limit(s=None):
    """Atom types the AutoDock-GPU binary of this job supports."""
    import anchored_docking as ad

    if s is not None and s.adg_path:
        for name in ad.ADGPU_BINARY_NAMES:
            candidate = os.path.join(s.adg_path, name)
            if os.path.isfile(candidate):
                return ad.atom_type_limit(candidate)
    return MAX_ATOM_TYPES


def disk_per_ligand_mb(size, anchors=3):
    """Intermediate files for one ligand: every anchor adds one grid map."""
    return (MB_PER_LIGAND_20A + MB_PER_ANCHOR_20A * (anchors - 3)) * (size / 20.0) ** 3


def anchors_in_use(s):
    """The most anchors any ligand will get."""
    if s.anchor_mode == "spread":
        most = max((r.anchors for r in s.kept if r.anchors), default=s.anchor_num)
        return min(s.anchor_num, most)
    return 3 if s.anchor_mode == "area" else s.anchor_num


def estimate_disk_mb(n_ligands, size, clean, anchors=3):
    per_ligand = disk_per_ligand_mb(size, anchors)
    if clean:
        return per_ligand + n_ligands * MB_KEPT_PER_LIGAND  # peak for one ligand plus what is kept
    return per_ligand * n_ligands


def nearest_existing_parent(path):
    path = os.path.abspath(path)
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:
            break
        path = parent
    return path


def validate_outdir(out, protected):
    """Return a message if `out` is unsafe or unusable as an output folder, otherwise None."""
    if re.search(r"\s", out):
        return "Please avoid spaces in the output folder path."
    real = os.path.realpath(out)
    prefix = real.rstrip(os.sep) + os.sep
    for p in protected:
        rp = os.path.realpath(p)
        if rp == real or rp.startswith(prefix):
            return f"That folder contains {p}, and an existing output folder is emptied. Choose another folder."
    if os.path.isfile(out):
        return "That path is a file, not a folder."
    parent = nearest_existing_parent(os.path.dirname(out))
    if not os.path.isdir(parent) or not os.access(parent, os.W_OK | os.X_OK):
        return f"You do not have permission to write in {parent}."
    if os.path.isdir(out) and not os.access(out, os.W_OK | os.X_OK):
        return f"You do not have permission to write in {out}."
    return None


def adfr_ready(path=None):
    """True if ADFRsuite's prepare_receptor is installed (the pixi 'adfr' environment)."""
    if path is None:
        import anchored_docking as ad

        path = ad.DEFAULT_ADFR_PATH
    exe = os.path.join(path, "bin", "prepare_receptor")
    return os.path.isfile(exe) and os.access(exe, os.X_OK)


def manual_core_anchors(smarts, template_mol):
    """Validate a manual core pattern the way anchored_docking.py will use it.

    Returns the 0-based atom positions that can be used as anchors; raises ValueError otherwise.
    """
    from rdkit import Chem, RDLogger

    import anchored_docking as ad

    RDLogger.DisableLog("rdApp.*")
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError("That is not a valid SMARTS pattern.")
    if Chem.MolFromSmiles(smarts) is None:
        raise ValueError("anchored_docking.py also reads this pattern as SMILES, so it must be valid SMILES "
                         "(no SMARTS-only wildcards or queries).")
    if not template_mol.HasSubstructMatch(query):
        raise ValueError("The pattern does not match the template ligand.")
    anchors = [idx for _, idx in ad.find_anchors(smarts)]
    if not anchors:
        raise ValueError("No unique anchor atoms were found in that pattern.")
    return anchors


def parse_anchor_indices(text, available):
    try:
        values = [int(t) for t in text.split()]
    except ValueError:
        raise ValueError("Enter whole numbers separated by spaces.")
    bad = [v for v in values if v not in available]
    if bad:
        raise ValueError(f"Not usable as anchors: {bad}. Usable positions: {sorted(available)}")
    return values


# --------------------------------------------------------------------------
# Input checks
# --------------------------------------------------------------------------


def analyse_receptor(path):
    info = {"atoms": 0, "hydrogens": 0, "waters": 0, "hetero": 0, "hetnames": set(),
            "models": 0, "altlocs": 0, "chains": set(), "coords": []}
    with open(path, errors="replace") as fh:
        for line in fh:
            record = line[:6]
            if record.startswith("MODEL"):
                info["models"] += 1
                continue
            if record not in ("ATOM  ", "HETATM"):
                continue
            resname = line[17:20].strip()
            name = line[12:16].strip()
            element = line[76:78].strip().upper() or re.sub(r"[^A-Za-z]", "", name)[:1].upper()
            if record == "HETATM" and resname in WATER_NAMES:
                info["waters"] += 1
                continue
            if line[16:17] not in (" ", "A", ""):
                info["altlocs"] += 1
                continue  # count each alternate location only once
            if record == "HETATM":
                info["hetero"] += 1
                info["hetnames"].add(resname)
            else:
                info["atoms"] += 1
                info["chains"].add(line[21:22])
            if element in ("H", "D"):
                info["hydrogens"] += 1
                continue
            try:
                info["coords"].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return info


def report_receptor(info, hydrogens_fixable=False, check_hetero=True):
    """Print findings. Returns False if the file is unusable.

    hydrogens_fixable: the wizard offers to add missing hydrogens, so only say they are missing.
    check_hetero: False for a receptor the wizard cleaned itself, where leftover molecules were chosen.
    """
    if info["atoms"] == 0:
        fail("No protein atoms (ATOM records) were found. Is this a PDB file?")
        return False
    chains = ", ".join(sorted(c.strip() or "-" for c in info["chains"]))
    ok(f"Protein with {info['atoms'] - info['hydrogens']} heavy atoms (chains: {chains}).")
    if info["hydrogens"] == 0 and hydrogens_fixable:
        warn("The receptor has no hydrogen atoms.", record=False)
    elif info["hydrogens"] == 0:
        warn("The receptor has NO hydrogen atoms. They are not added automatically, so the "
             "results will be poor. Add them first (e.g. `h_add` in PyMOL) and save a new PDB.")
    else:
        ok(f"Hydrogens present ({info['hydrogens']}).")
    if info["hetero"] and check_hetero:
        names = ", ".join(sorted(info["hetnames"]))
        warn(f"The file also contains non-protein molecules ({names}). Waters are removed "
             "automatically, but leftovers such as the original ligand can interfere. "
             "Remove them unless you need them.")
    if info["models"] > 1:
        warn(f"The file contains {info['models']} models (e.g. an NMR ensemble). Keep just one.")
    if info["altlocs"]:
        warn(f"{info['altlocs']} atoms have alternate locations. Keep one conformation.")
    return True


def load_template(path):
    from rdkit import Chem

    try:
        mol = Chem.SDMolSupplier(path)[0]
    except (IndexError, OSError):
        mol = None
    return mol


def count_molecules(path):
    """Number of records in an SDF file (counts '$$$$' terminators)."""
    with open(path, errors="replace") as fh:
        return sum(1 for line in fh if line.startswith("$$$$"))


def report_template(mol, receptor_coords, n_molecules=1):
    """Print findings. Returns the template's longest dimension in Angstrom, or None if unusable."""
    import numpy as np
    from scipy.spatial import cKDTree

    if mol is None:
        fail("RDKit could not read a molecule from this SDF file.")
        return None
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        fail("The template has no 3D coordinates. It must be a ligand positioned in the receptor.")
        return None
    positions = mol.GetConformer().GetPositions()
    if np.ptp(positions, axis=0).min() < 0.01:
        fail("All template atoms lie in one plane, so these are 2D coordinates. "
             "The template must be a 3D ligand positioned in the receptor.")
        return None
    ok(f"Template: {mol.GetNumHeavyAtoms()} heavy atoms.")
    if n_molecules > 1:
        warn(f"The template file contains {n_molecules} molecules; only the first is used. "
             "Put the one you want first, or use a single-molecule file.")
    heavy = np.array([p for p, a in zip(positions, mol.GetAtoms()) if a.GetAtomicNum() > 1])
    if len(receptor_coords):
        distances, _ = cKDTree(np.array(receptor_coords)).query(heavy)
        close = float(np.mean(distances <= 4.5))
        if distances.min() < 1.0:
            warn("Template atoms sit directly on top of receptor atoms. Is the ligand still "
                 "inside the receptor PDB file, or is it clashing with the protein?")
        elif close < 0.5:
            warn(f"Only {close:.0%} of the template atoms are near the receptor. The template and "
                 "receptor are probably in different coordinate systems. Export both from the "
                 "same structure without moving either one.")
        else:
            ok(f"Template sits in the receptor ({close:.0%} of its atoms are within 4.5 A of protein atoms).")
    extent = float((heavy.max(axis=0) - heavy.min(axis=0)).max())
    note(f"The template spans about {extent:.1f} A at its longest.")
    return extent


@dataclass
class LigandResult:
    index: int
    name: str
    mol: object
    usable: bool
    reason: str = ""
    anchors: int = 0
    core_atoms: int = 0
    notes: list = field(default_factory=list)  # problems that do not make the ligand impossible, but likely
    atom_types: int = 0     # AutoDock atom types needed with every usable core atom as an anchor (0 = not measured)


def load_atom_types():
    import json

    with open(os.path.join(REPO_DIR, "scripts", "base_parameters.json")) as fh:
        return json.load(fh)["ad4"]


def untyped_symbols(mol, atom_types):
    """Element symbols of heavy atoms that no AutoDock atom type pattern covers."""
    from rdkit import Chem

    typed = set()
    for info in atom_types:
        pattern = Chem.MolFromSmarts(info["smarts"])
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern):
            typed.update(match)
    return sorted({a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetIdx() not in typed})


def ligand_notes(mol, core_atoms, atom_types):
    """Likely problems for a ligand that has a shared core with the template."""
    from rdkit import Chem

    notes = []
    if core_atoms < MIN_CORE_ATOMS:
        notes.append(f"shares only {core_atoms} heavy atom(s) with the template; it is probably from a "
                     "different chemical series, so the anchors may not make sense")
    fragments = len(Chem.GetMolFrags(mol))
    if fragments > 1:
        notes.append(f"has {fragments} separate fragments (a salt or counter-ion?); remove the extra "
                     "fragments before docking")
    symbols = untyped_symbols(mol, atom_types)
    if symbols:
        notes.append(f"contains atoms AutoDock has no atom type for ({', '.join(symbols)}); docking fails "
                     "if one of them is chosen as an anchor, and may fail earlier")
    return notes


def ligand_name(mol, index):
    if mol is not None and mol.HasProp("_Name") and mol.GetProp("_Name").strip():
        return mol.GetProp("_Name").strip()
    return f"molecule {index + 1}"


def check_ligands(template_mol, path):
    """Run the same core/anchor search the docking script uses, for every ligand."""
    from rdkit import Chem, RDLogger

    import anchored_docking as ad

    RDLogger.DisableLog("rdApp.*")
    mols = list(Chem.SDMolSupplier(path))
    atom_types = load_atom_types()
    results = []
    total = len(mols)
    try:
        for i, mol in enumerate(mols):
            name = ligand_name(mol, i)
            sys.stdout.write(f"\r  Checking ligand {i + 1} of {total}...   ")
            sys.stdout.flush()
            if mol is None:
                results.append(LigandResult(i, name, None, False, "RDKit could not read this molecule"))
                continue
            try:
                search = ad.find_anchor_points(template_mol, mol, 0, 100, None)
            except Exception as exc:  # same code path as the real run; report instead of crash
                results.append(LigandResult(i, name, mol, False, f"could not find a shared core ({type(exc).__name__})"))
                continue
            if not search.anchor_patterns:
                results.append(LigandResult(i, name, mol, False, "no shared core with the template"))
                continue
            try:
                ad.embed_ligand(mol)  # the same 3D generation the job will do
            except ad.LigandError as exc:
                results.append(LigandResult(i, name, mol, False, str(exc)))
                continue
            core_atoms = search.pattern_mol.GetNumHeavyAtoms()
            results.append(LigandResult(i, name, mol, True, anchors=len(search.anchor_patterns),
                                        core_atoms=core_atoms,
                                        notes=ligand_notes(mol, core_atoms, atom_types)))
    except KeyboardInterrupt:
        say("\n  Check interrupted: the remaining ligands were not checked.")
        for i in range(len(results), total):
            results.append(LigandResult(i, ligand_name(mols[i], i), mols[i], mols[i] is not None,
                                        "" if mols[i] is not None else "RDKit could not read this molecule",
                                        anchors=0, core_atoms=0))
    sys.stdout.write("\r" + " " * 40 + "\r")
    return results


def flagged_ligands(results):
    """Ligands that cannot be docked, or are likely to cause trouble."""
    return [r for r in results if not r.usable or r.notes]


def report_ligands(results):
    good = [r for r in results if r.usable]
    bad = [r for r in results if not r.usable]
    if not results:
        fail("No molecules were found in this file.")
        return
    ok(f"{len(results)} ligand(s) in the file; {len(good)} share a core with the template.")
    checked = [r for r in good if r.anchors]
    if checked:
        note(f"Smallest shared core: {min(r.core_atoms for r in checked)} heavy atoms. "
             f"Fewest usable anchor atoms in one ligand: {min(r.anchors for r in checked)}.")
    for r in bad[:15]:
        fail(f"#{r.index + 1} {r.name[:40]}: {r.reason}")
    if len(bad) > 15:
        note(f"... and {len(bad) - 15} more.")
    noted = [r for r in good if r.notes]
    for r in noted[:10]:
        for text in r.notes:
            warn(f"#{r.index + 1} {r.name[:40]} {text}.", record=False)
    if len(noted) > 10:
        warn(f"{len(noted) - 10} more ligands have similar problems.", record=False)


# --------------------------------------------------------------------------
# Settings, command and job script
# --------------------------------------------------------------------------


@dataclass
class Settings:
    receptor: str = ""
    template: str = ""
    ligands: str = ""
    kept: list = field(default_factory=list)  # LigandResult objects that will be docked
    all_count: int = 0
    outdir: str = ""
    force: bool = False
    anchor_num: int = 3
    anchor_mode: str = "random"
    size: float = 20.0
    samples: int = 100
    seed: object = None
    query_smarts: object = None
    anchor_indices: object = None
    clean: bool = False
    verbose: bool = False
    mode: str = "submit"
    jobname: str = "anchored_dock"
    partition: str = DEFAULT_PARTITION
    reservation: str = ""
    gpus: int = 1
    time_limit: str = "0:30:00"
    modules: str = DEFAULT_MODULES
    pixi_dir: str = ""
    ligands_used: str = ""
    add_hydrogens: bool = False
    adg_path: str = ""          # folder of a custom AutoDock-GPU build with a higher atom-type limit
    protonated: object = None   # protein_prep.ProteinPrep when the wizard protonated the receptor
    original_receptor: str = ""  # the user's own receptor file, before protonation (separate-files route)
    template_view: str = ""     # the template with hydrogens, for display only (complex route)
    prepared: object = None     # structure_prep.Prepared when the receptor and template come from a complex
    structure_path: str = ""
    resume: bool = False
    remaining: int = 0          # ligands the job will actually dock (fewer than len(kept) when resuming)
    script_path: str = ""
    template_mol: object = None
    extent: float = 0.0
    anchor_limit: int = 0


def jobdir(s):
    return os.path.dirname(s.outdir)


def build_args(s):
    args = ["--template_file", s.template, "--input_file", s.ligands_used or s.ligands,
            "--receptor_file", s.receptor, "--docking_dir", s.outdir,
            "--scripts_path", os.path.join(REPO_DIR, "scripts"),
            "--anchor_num", str(s.anchor_num), "--anchor_mode", s.anchor_mode,
            "--size", str(s.size), "--samples", str(s.samples)]
    if s.adg_path:
        args += ["--adg_path", s.adg_path]
    if s.query_smarts:
        args += ["--query_smarts", s.query_smarts]
    if s.anchor_indices:
        args += ["--anchor_indices"] + [str(i) for i in s.anchor_indices]
    if s.seed is not None:
        args += ["--random_seed", str(s.seed)]
    if s.prepared is not None:
        # A template taken from a structure file has no bond orders, so RDKit must not guess its hydrogens.
        args += ["--template_view_file", s.template_view or prepared_paths(s)[2]]
    for flag, on in (("--force", s.force and not s.resume), ("--resume", s.resume),
                     ("--add_hydrogens", s.add_hydrogens), ("--clean", s.clean), ("--verbose", s.verbose)):
        if on:
            args.append(flag)
    return args


def multiline_command(prefix, args):
    groups, current = [], []
    for token in args:
        if token.startswith("--") and current:
            groups.append(current)
            current = []
        current.append(token)
    groups.append(current)
    lines = [prefix] + ["    " + " ".join(shlex.quote(t) for t in g) for g in groups]
    return " \\\n".join(lines)


MAIL_DOMAIN = "crick.ac.uk"


def mail_address():
    """Where SLURM sends the notification: the current user name plus the lab's mail domain (None if unknown)."""
    try:
        user = getpass.getuser()
    except Exception:
        return None
    return f"{user}@{MAIL_DOMAIN}" if re.fullmatch(r"[A-Za-z0-9._-]+", user) else None


def build_script(s):
    mail = mail_address()
    lines = ["#!/bin/bash",
             f"#SBATCH --partition={s.partition}",
             *([f"#SBATCH --reservation={s.reservation}"] if s.reservation else []),
             f"#SBATCH --gres=gpu:{s.gpus}",
             f"#SBATCH --job-name={s.jobname}",
             f"#SBATCH --output={os.path.join(jobdir(s), s.jobname)}_%j.log",
             f"#SBATCH --time={s.time_limit}",
             *(["#SBATCH --mail-type=END,FAIL", f"#SBATCH --mail-user={mail}"] if mail else []),
             ""]
    for module in s.modules.split(","):
        if module.strip():
            lines.append(f"module load {module.strip()}")
    lines += ["", f'export PATH="{s.pixi_dir}:$PATH"', f"cd {shlex.quote(REPO_DIR)}", "",
              multiline_command("pixi run python anchored_docking.py", build_args(s)), ""]
    return "\n".join(lines)


def write_ligands(path, kept):
    from rdkit import Chem

    with Chem.SDWriter(path) as writer:
        for r in kept:
            writer.write(r.mol)


def prepare_ligand_file(s):
    """If some ligands were excluded, write a reduced SDF for the job to use."""
    if len(s.kept) == s.all_count:
        s.ligands_used = s.ligands
        return
    s.ligands_used = os.path.join(jobdir(s), f"{s.jobname}_ligands.sdf")
    write_ligands(s.ligands_used, s.kept)


def prepared_paths(s):
    """Where the cleaned receptor, the extracted template and its display version are saved for the job."""
    return (os.path.join(jobdir(s), f"{s.jobname}_receptor.pdb"),
            os.path.join(jobdir(s), f"{s.jobname}_template.sdf"),
            os.path.join(jobdir(s), f"{s.jobname}_template_with_hydrogens.sdf"))


def finalize_structure_files(s):
    """Save the cleaned receptor and template next to the results and point the job at them.

    Until now they lived in a temporary folder that disappears when the assistant exits; the job
    (and a later resume) needs files that stay.
    """
    if s.prepared is None and s.protonated is None:
        return
    os.makedirs(jobdir(s), exist_ok=True)
    receptor, template, template_view = prepared_paths(s)
    if s.prepared is not None:
        with open(receptor, "w") as fh:
            fh.write(s.prepared.receptor_pdb)       # already the protonated text if the wizard protonated it
        with open(template, "w") as fh:
            fh.write(s.prepared.template_molblock)
        with open(template_view, "w") as fh:  # without a dictionary entry this is the heavy-atom template itself
            fh.write(s.prepared.template_view_molblock or s.prepared.template_molblock)
        s.receptor, s.template, s.template_view = receptor, template, template_view
    else:  # separate files: only the receptor is new, the user's template stays where it is
        with open(receptor, "w") as fh:
            fh.write(s.protonated.pdb)
        s.receptor = receptor
    if s.protonated is not None:
        with open(os.path.join(jobdir(s), f"{s.jobname}_protonation_report.txt"), "w") as fh:
            fh.write(s.protonated.report_text())


def previous_run(s, out):
    """Compare a folder's earlier run with the current choices.

    Returns None if the folder holds no run of this workflow, otherwise
    (names of settings that differ, number of ligands already docked).
    """
    import tempfile

    import anchored_docking as ad

    previous = ad.read_run_settings(out)
    if previous is None:
        return None
    ligands, temp = s.ligands, None
    if len(s.kept) != s.all_count:  # the job will use the reduced list, so compare against that
        handle, temp = tempfile.mkstemp(suffix=".sdf")
        os.close(handle)
        write_ligands(temp, s.kept)
        ligands = temp
    try:
        current = ad.run_signature(s.template, ligands, s.receptor, s.anchor_num, s.anchor_mode, s.size,
                                   s.samples, s.seed, s.query_smarts, s.anchor_indices, s.add_hydrogens)
    finally:
        if temp:
            os.remove(temp)
    done = sum(os.path.isfile(os.path.join(out, str(i), "ligand_docked.sdf")) for i in range(len(s.kept)))
    return ad.signature_differences(previous, current), done


# --------------------------------------------------------------------------
# Wizard steps
# --------------------------------------------------------------------------


def offer_hydrogens(s, info):
    """If the receptor has no hydrogens, offer to have ADFR add them."""
    if info["hydrogens"] > 0:
        return
    say("\n  Docking needs hydrogens on the receptor. ADFR can add them with a quick rule-based guess")
    say("  (histidine protonation is not optimised). Receptors that already have hydrogens are not touched.")
    s.add_hydrogens = ask_yes_no("Add hydrogens automatically?", True)
    if s.add_hydrogens:
        warn("Hydrogens will be added to the receptor by a quick guess; check histidine protonation if it matters "
             "for your target.")
    else:
        warn("The receptor has no hydrogens and none will be added; the results will be poor.")


def stage_receptor(s, text):
    """Make `text` the receptor the rest of the wizard looks at (in a temporary folder, saved with the job later)."""
    import atexit
    import tempfile

    folder = tempfile.mkdtemp(prefix="anchored_docking_receptor_")
    atexit.register(shutil.rmtree, folder, ignore_errors=True)
    s.receptor = os.path.join(folder, "receptor.pdb")
    with open(s.receptor, "w") as fh:
        fh.write(text)


def step_protein(s, info, template_mol):
    """Protein preparation: report gaps and missing atoms, then protonate at a chosen pH if there are no hydrogens."""
    import numpy as np

    import protein_prep as pp

    say("\n  Protein preparation")
    with open(s.receptor, errors="replace") as fh:
        text = fh.read()
    template_xyz = np.array([p for p, a in zip(template_mol.GetConformer().GetPositions(), template_mol.GetAtoms())
                             if a.GetAtomicNum() > 1])
    notes, urgent = pp.quality_lines(pp.read_residues(text), template_xyz)
    for line in notes:
        note(line)
    for line in urgent:
        warn(line)
    if not notes and not urgent:
        ok("No gaps in the chains and no residues with missing heavy atoms.")

    if info["hydrogens"] > 0:
        ok(f"The receptor already has hydrogens ({info['hydrogens']}); they are used as they are.")
        note("Protonation states (histidine form, charged or neutral acids) are not checked. To have them "
             "predicted, start from a receptor without hydrogens.")
        return

    say("\n  Docking needs hydrogens on the receptor. PDB2PQR with PROPKA predicts the charge state of each")
    say("  residue at a pH you choose (histidine form, aspartate and glutamate, termini), places hydrogens so")
    say("  that hydrogen bonds are satisfied, and rebuilds missing side-chain atoms. It takes about a minute.")
    if ask_yes_no("Protonate the protein with PDB2PQR and PROPKA?", True):
        ph = ask_number("pH", pp.DEFAULT_PH, minimum=0.0, maximum=14.0)
        say("  Working...")
        try:
            prep = pp.protonate(text, template_xyz, ph)
        except pp.ProteinPrepError as exc:
            warn(f"Protonation with PDB2PQR did not work: {exc}")
            note("Falling back to ADFR's quicker rule-based hydrogen placement.")
            offer_hydrogens(s, info)
            return
        s.protonated = prep
        s.original_receptor = s.receptor
        stage_receptor(s, prep.pdb)
        if s.prepared is not None:
            s.prepared.receptor_pdb = prep.pdb
        for line in prep.lines():
            ok(line)
        notable = [p for p in prep.pocket if p.notable]
        if notable:
            say(f"\n  Residues within {pp.POCKET_DISTANCE:g} A of the template that deserve a look:")
            table([(p.label, f"{p.pka:.1f}", p.state) for p in notable], headers=("Residue", "pKa", "State"))
        elif prep.pocket:
            note(f"The {len(prep.pocket)} titratable residue(s) near the template all keep their usual state.")
        for problem in prep.problems:
            warn(f"PDB2PQR: {problem}.")
        if prep.kept_without_hydrogens:
            warn(f"{', '.join(prep.kept_without_hydrogens)} kept in the receptor got no hydrogens (PDB2PQR only "
                 "protonates protein).")
        note("The full record is saved with the results as <job name>_protonation_report.txt.")
    else:
        offer_hydrogens(s, info)


def show_inventory(inv):
    chains = ", ".join(f"{name} ({n} residues)" for name, n in inv.protein_chains.items())
    ok(f"Protein chain(s): {chains}.")
    if inv.nucleic_chains:
        note(f"Nucleic acid chain(s) {', '.join(inv.nucleic_chains)} will not be used.")
    others = [h for h in inv.hets if h.kind != "sugar"]
    kinds = {}
    for h in others:
        kinds.setdefault(h.kind, []).append(h.name)
    parts = [f"{len(names)} {kind}{'s' if len(names) > 1 and kind != 'additive' else ''} ({', '.join(sorted(set(names)))})"
             for kind, names in kinds.items()]
    note(f"Also in the file: {inv.waters} water(s)" + (f", {', '.join(parts)}" if parts else "")
         + (f", {sum(1 for h in inv.hets if h.kind == 'sugar')} sugar residue(s)" if any(h.kind == "sugar" for h in inv.hets) else "")
         + ".")


def molecule_rows(shown, components):
    """Table rows for the molecules offered as template: what each one is, in terms a chemist recognises.

    The name and complete formula come from the PDB's chemical component dictionary. Without them (offline, or an
    unregistered molecule) the formula of the heavy atoms in the file is shown instead, marked as having no hydrogens.
    """
    rows = []
    for n, het in enumerate(shown, 1):
        info = components.describe(het.name)
        formula = info.formula or (f"{het.heavy_formula} (no H)" if het.heavy_formula else "")
        rows.append((n, info.name or "(name unavailable)", formula, het.name, het.chain.strip() or "-", het.number,
                     het.kind, het.heavy, het.contacts))
    return rows


def per_molecule_notes(shown, components, pick):
    """One line per distinct molecule instead of one per row: copies of a molecule (the same residue code in
    several chains) share their name and SMILES, so they are listed once under all their row numbers.

    pick(info, het) returns the text for a molecule, or "" for none. Returns [("1, 2", text), ...] in table order.
    """
    grouped = {}
    for n, het in enumerate(shown, 1):
        text = pick(components.describe(het.name), het)
        if text:
            grouped.setdefault((het.name, text), []).append(str(n))
    return [(", ".join(numbers), text) for (_, text), numbers in grouped.items()]


def choose_template(inv, components=None):
    """Ask which molecule is the template. Returns a structure_prep.Het, or None to give up on this file."""
    import structure_prep as sp

    components = components or sp.Components()

    candidates = sp.template_candidates(inv)
    if not candidates:
        fail("No molecule in this file is big enough to be a ligand, so there is no template to extract.")
        return None
    suggestion = sp.suggested_template(inv)
    say("\n  Which molecule is the template ligand (the one to anchor to)?")
    shown = candidates[:12]
    found = components.prefetch([h.name for h in shown])
    table(molecule_rows(shown, components),
          headers=("#", "Molecule", "Formula", "Code", "Chain", "Residue No.", "Kind", "Heavy atoms",
                    "Protein contacts"))
    if shown and not found:
        note("The names could not be looked up at the PDB (no connection?), so only the formulas of the atoms in the "
             "file are shown. Hydrogens are not in the file, so the real formulas have more H.")
    full = per_molecule_notes(shown, components, lambda info, het: info.systematic
                              if info.systematic and info.systematic.lower() != info.name.lower() else "")
    if full:
        say("\n  Full chemical names:")
        for numbers, text in full:
            say(f"    {numbers}  {text}")
    smiles = per_molecule_notes(shown, components,
                                lambda info, het: info.smiles if het.kind in ("ligand", "cofactor") else "")
    if smiles:
        say("\n  SMILES of the ligand-like molecules (paste one into a drawing tool such as PubChem Sketcher to see it):")
        for numbers, text in smiles:
            say(f"    {numbers}  {text}")
    if len(candidates) > len(shown):
        note(f"{len(candidates) - len(shown)} smaller or less likely molecules are not shown.")
    default = str(shown.index(suggestion) + 1) if suggestion in shown else None
    if default is None:
        warn("None of the molecules looks like a drug-like ligand; pick the one you want, if any.", record=False)

    def pick(text):
        if text.lower() in ("q", "none"):
            return None
        if text.isdigit() and 1 <= int(text) <= len(shown):
            return shown[int(text) - 1]
        raise ValueError(f"Enter a number from 1 to {len(shown)}, or q to choose another file.")

    return ask("Template number (q = choose another file)", default, pick)


def choose_chains(inv, distances):
    """Ask which protein chains to keep (only when there is a choice)."""
    import structure_prep as sp

    defaults = sp.default_chains(distances)
    if len(inv.protein_chains) == 1:
        return list(inv.protein_chains)
    say("\n  The structure has several protein chains. Chains touching the template's neighbourhood are")
    say(f"  suggested (an atom within {sp.CHAIN_NEAR:g} A of the template).")
    table([(name, n, f"{distances[name]:.1f}", "keep" if name in defaults else "") for name, n in inv.protein_chains.items()],
          headers=("Chain", "Residues", "Nearest to template (A)", "Suggested"))

    def validate(text):
        names = text.replace(",", " ").split()
        unknown = [n for n in names if n not in inv.protein_chains]
        if unknown or not names:
            raise ValueError(f"Choose from: {', '.join(inv.protein_chains)} (separated by spaces).")
        return sorted(set(names), key=list(inv.protein_chains).index)

    return ask("Chains to keep", " ".join(defaults), validate)


def ask_molecules(plan):
    """Ask, for each ion, cofactor or other molecule near the template, whether to keep it in the receptor."""
    import structure_prep as sp

    if not plan.ask:
        return []
    say("\n  Waters, crystallisation additives and sugars are removed. For the molecules below, decide whether")
    say("  they belong to the receptor (for example a catalytic metal or a cofactor in the pocket).")
    kept = []
    if len(plan.ask) > sp.GROUP_ABOVE:
        groups = {}
        for het in plan.ask:
            groups.setdefault(het.name, []).append(het)
        say(f"  There are {len(plan.ask)} of them, so molecules of the same kind are asked about together.")
        for name, members in groups.items():
            nearest = min(h.distance for h in members)
            if ask_yes_no(f"Keep all {len(members)} {name} ({members[0].kind}, nearest {nearest:.1f} A from the template)?",
                          nearest <= sp.KEEP_DEFAULT_RADIUS):
                kept.extend(members)
        return kept
    for het in plan.ask:
        text = f"Keep {het.label} ({het.kind}, {het.heavy} atom{'s' if het.heavy != 1 else ''}, {het.distance:.1f} A from the template)?"
        if ask_yes_no(text, sp.default_keep(het)):
            kept.append(het)
    return kept


def stage_files(s, prepared):
    """Write the cleaned receptor and template to a temporary folder so the rest of the wizard can use them."""
    import atexit
    import tempfile

    folder = tempfile.mkdtemp(prefix="anchored_docking_")
    atexit.register(shutil.rmtree, folder, ignore_errors=True)
    s.receptor = os.path.join(folder, "receptor.pdb")
    s.template = os.path.join(folder, "template.sdf")
    with open(s.receptor, "w") as fh:
        fh.write(prepared.receptor_pdb)
    with open(s.template, "w") as fh:
        fh.write(prepared.template_molblock)
    s.prepared = prepared


def prepare_from_complex(s):
    """Guided clean-up of a protein-ligand complex. Sets s.receptor, s.template and s.prepared."""
    import structure_prep as sp

    say("  A PDB (.pdb) or mmCIF (.cif) file, for example downloaded from rcsb.org. It is only read: your")
    say("  file is not changed. Press Tab to complete file names.\n")
    components = sp.Components()   # what the molecules are, looked up once and shared by the steps below
    while True:
        s.structure_path = ask("Structure file", validate=file_validator(sp.STRUCTURE_EXTENSIONS),
                               completer=path_completer)
        try:
            structure, n_models, altlocs = sp.load_structure(s.structure_path)
            inv = sp.inventory(structure, n_models, altlocs, s.structure_path)
        except sp.PrepError as exc:
            fail(str(exc))
            continue
        say("")
        show_inventory(inv)
        template = choose_template(inv, components)
        if template is None:
            continue
        break
    distances = sp.measure_distances(inv, template)
    chains = choose_chains(inv, distances)
    plan = sp.plan_molecules(inv, template, chains)
    keep = ask_molecules(plan)
    say(f"\n  Looking up the chemistry of {template.name} at the PDB to give the template its hydrogens for display...")
    try:
        prepared = sp.prepare(structure, inv, template, chains, keep, ccd=components.fetch)
    except sp.PrepError as exc:
        fail(str(exc))
        raise SystemExit(1)
    say("\n  What was done with the structure:")
    for line in prepared.summary:
        say(f"    - {line}")
    if keep:
        warn(f"Kept in the receptor: {', '.join(h.label for h in keep)}. Ions and cofactors may lack AutoDock "
             "partial charges (ADFR reports this for metals such as zinc), which weakens their effect.")
    stage_files(s, prepared)
    return inv, template


def step_files(s):
    heading(1, "Protein structure")
    say("  What do you have?")
    say("    1. A structure file of the protein-ligand complex (PDB or mmCIF). The assistant cleans it up:")
    say("       it extracts the template ligand and removes waters, additives and other chains.")
    say("    2. A cleaned protein PDB file and a separate template ligand SDF file.")

    def pick(text):
        if text in ("1", "2"):
            return text
        raise ValueError("Please enter 1 or 2.")

    from_complex = ask("Choose a number", "1", pick) == "1"
    say("")
    mark = len(WARNINGS)
    if from_complex:
        prepare_from_complex(s)
        info = analyse_receptor(s.receptor)
        say("")
        report_receptor(info, hydrogens_fixable=True, check_hetero=False)
    else:
        say("  A PDB file of the protein, without the ligand. Press Tab to complete file names.\n")
        while True:
            del WARNINGS[mark:]  # forget warnings about a file the user is replacing
            s.receptor = ask("Receptor PDB file", validate=file_validator((".pdb",)), completer=path_completer)
            info = analyse_receptor(s.receptor)
            if report_receptor(info, hydrogens_fixable=True):
                break
    receptor_coords = info["coords"]

    heading(2, "Template ligand")
    mark = len(WARNINGS)
    if from_complex:
        say(f"  Taken from the structure: {s.prepared.template.label}.\n")
        template_mol = load_template(s.template)
        extent = report_template(template_mol, receptor_coords, 1)
        if extent is None:
            fail("The extracted template could not be used. Choose the separate-files option instead.")
            raise SystemExit(1)
    else:
        say("  A ligand that is already in the right place in the receptor, e.g. from a crystal")
        say("  structure (SDF file). Only the first molecule in the file is used.\n")
        while True:
            del WARNINGS[mark:]
            s.template = ask("Template ligand SDF file", validate=file_validator((".sdf",)), completer=path_completer)
            template_mol = load_template(s.template)
            extent = report_template(template_mol, receptor_coords, count_molecules(s.template))
            if extent is not None:
                break
    s.template_mol = template_mol
    s.extent = extent
    step_protein(s, info, template_mol)

    heading(3, "Ligands to dock")
    say("  An SDF file with one or more ligands that share a core with the template.")
    say("  Each ligand is docked separately.\n")
    while True:
        s.ligands = ask("Ligand SDF file", validate=file_validator((".sdf",)), completer=path_completer)
        say("")
        results = check_ligands(template_mol, s.ligands)
        report_ligands(results)
        good = [r for r in results if r.usable]
        if not good:
            fail("None of these ligands can be docked with this template. Choose another file.")
            continue
        break
    s.all_count = len(results)
    s.kept = good
    flagged = flagged_ligands(results)
    if flagged:
        flagged_ids = {id(r) for r in flagged}
        remaining = [r for r in results if id(r) not in flagged_ids]
        say("")
        if remaining and ask_yes_no(f"Leave out the {len(flagged)} flagged ligand(s) and dock the other {len(remaining)}?", True):
            s.kept = remaining
            say(f"  {len(remaining)} ligand(s) will be docked.")
        else:
            s.kept = results
            warn("Flagged ligands are included. Ones that cannot be docked are reported as failed; "
                 "the other ligands are still docked.")
    s.remaining = len(s.kept)
    s.anchor_limit = min((r.anchors for r in s.kept if r.anchors), default=0)


def standard_ligand_types():
    """Number of ligand atom types AutoGrid makes maps for (read from scripts/write-gpf.py)."""
    try:
        with open(os.path.join(REPO_DIR, "scripts", "write-gpf.py"), errors="replace") as fh:
            match = re.search(r"^ligand_types (.+)$", fh.read(), re.M)
        return len(match.group(1).split())
    except (OSError, AttributeError):
        return 15


def atom_types_needed(template_mol, mol):
    """Exactly how many AutoDock atom types the ligand needs with every usable core atom as an anchor.

    This runs the same anchor choice and Meeko preparation as the job, in a temporary folder, and counts the
    types in the PDBQT. Returns None if the ligand cannot be prepared (it is reported elsewhere).
    """
    import argparse
    import tempfile

    import anchored_docking as ad

    args = argparse.Namespace(anchor_num=0, samples=100, random_seed=None, query_smarts=None,
                              anchor_indices=None, anchor_mode="spread", verbose=False)
    tools = {"scripts_path": os.path.join(REPO_DIR, "scripts"), "atom_type_limit": 10 ** 9}
    try:
        with tempfile.TemporaryDirectory(prefix="atom_types_") as folder:
            ad.prepare_anchored_ligand(args, tools, template_mol, mol, folder)
            return len(ad.ligand_atom_types(os.path.join(folder, "ligand.pdbqt")))
    except Exception:  # noqa: BLE001 - a ligand that cannot be prepared fails later with its own message
        return None


def measure_atom_types(s):
    """Fill in LigandResult.atom_types for the ligands that could need more atom types than the default AutoDock-GPU has.

    Other atoms of a ligand have at most standard_ligand_types() different types, so a ligand whose usable core
    atoms number MAX_ATOM_TYPES minus that always fits and is not measured.
    """
    from concurrent.futures import ThreadPoolExecutor

    always_fits = MAX_ATOM_TYPES - standard_ligand_types()
    todo = [r for r in s.kept if r.mol is not None and r.anchors > always_fits]
    if not todo:
        return
    say(f"\n  Checking how many AutoDock atom types {len(todo)} ligand(s) with large cores need...")
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 1)) as pool:
        for ligand, needed in zip(todo, pool.map(lambda r: atom_types_needed(s.template_mol, r.mol), todo)):
            ligand.atom_types = needed or 0
    unmeasured = [r for r in todo if not r.atom_types]
    if unmeasured:
        warn(f"The atom types of {len(unmeasured)} of {len(todo)} ligand(s) with large cores could not be measured, so "
             "the AutoDock-GPU with a higher atom-type limit was not considered for them. A ligand that needs "
             "more atom types than the cluster's AutoDock-GPU supports then gets fewer anchors, and the results "
             "table says so.")


def choose_autodock_gpu(s):
    """Use the AutoDock-GPU with a higher atom-type limit when some ligand needs it, without asking.

    The default AutoDock-GPU supports MAX_ATOM_TYPES atom types per ligand; the anchors count towards that
    limit. Ligands over it would get fewer anchors, so a build with a higher limit is used if one exists.
    """
    over = [r for r in s.kept if r.atom_types > MAX_ATOM_TYPES]
    if not over:
        return
    most = max(r.atom_types for r in over)
    build = find_custom_build()
    if build is None:
        warn(f"{len(over)} ligand(s) need more AutoDock atom types than the cluster's AutoDock-GPU supports "
             f"({most} against {MAX_ATOM_TYPES}). They will get fewer anchors; the results table says how many. "
             "scripts/build_autodock_gpu.sh builds a version with a higher limit.")
        return
    folder, limit = build
    s.adg_path = folder
    say(f"\n  {len(over)} ligand(s) need more AutoDock atom types than the cluster's AutoDock-GPU supports "
        f"({most} against {MAX_ATOM_TYPES}),")
    say(f"  so the build with a limit of {limit} atom types is used: {folder}")
    have = build_targets(folder)
    if have:
        names = {"80": "A100 (ga100)", "90": "H100 (gh100)"}
        note("It was built for: " + ", ".join(names.get(t, f"compute capability {t}") for t in have) + ".")
    if limit < most:
        warn(f"Even that build supports only {limit} atom types, fewer than the {most} the largest ligand needs.")


def step_settings(s):
    heading(4, "Docking settings")
    usable = s.anchor_limit
    measure_atom_types(s)
    choose_autodock_gpu(s)
    limit = atom_type_limit(s)
    default_anchors = min(3, usable, limit - 1) if usable else 3   # for the random mode
    say("  Anchors are atoms of the shared core (the maximum common substructure of each ligand and the")
    say("  template) that are pulled towards the template position.")
    if usable and usable < 3:
        warn(f"Only {usable} usable anchor atom(s) in the weakest ligand: the core will be held loosely there.")
    if usable:
        counts = [r.anchors for r in s.kept if r.anchors]
        used = f"{min(counts)}" if min(counts) == max(counts) else f"{min(counts)} to {max(counts)}"
        say(f"\n  Recommended: every usable core atom as an anchor ({used} per ligand), chosen to be spread across the")
        say("  shared core, 20 A search box.")
    else:
        say("\n  Recommended: every usable core atom as an anchor, spread across the shared core, 20 A search box.")
    say(f"  AutoDock-GPU supports {limit} atom types in all (each anchor is one, and so are the types of the other")
    say("  atoms). A ligand that would need more gets fewer anchors, and the results table says so.")
    s.anchor_mode, s.anchor_num = "spread", limit
    if not ask_yes_no("Use the recommended docking settings?", True):
        say("  'spread' uses up to the number of anchors you allow, always the ones furthest apart from each other, so")
        say("  the whole core is held in place. 'random' picks a number of anchors at random. 'area' tries many")
        say("  random sets and keeps the three that are furthest apart (always exactly 3 anchors).")
        s.anchor_mode = ask_choice("Anchor selection", ["spread", "random", "area"], "spread")
        if s.anchor_mode == "area" and usable and usable < 3:
            warn("The 'area' mode needs at least 3 usable anchors; using 'random' instead.")
            s.anchor_mode = "random"
        if s.anchor_mode == "spread":
            say("\n  A ligand whose core has fewer usable atoms than the number you give uses all of them.")
            s.anchor_num = ask_number("Maximum number of anchor atoms", limit, int, 1, maximum=limit)
        elif s.anchor_mode == "area":
            s.anchor_num = 3
            say("\n  More random attempts find a better-spread set of anchors but take longer.")
            s.samples = ask_number("Number of random attempts", 100, int, 1)
        else:
            say(f"\n  Each anchor adds an atom type and AutoDock-GPU supports {limit}, so a ligand with many anchors")
            say("  and many other atom types can exceed it; that ligand then fails with an explanation.")
            s.anchor_num = ask_number("Number of anchor atoms", default_anchors, int, 1, maximum=limit - 1)
            if usable and s.anchor_num > usable:
                warn(f"The weakest ligand has only {usable} usable anchor atom(s); the job will fail on it.")
        say("\n  The search box is the region around the template where the ligand may move.")
        s.size = ask_number("Search box size in Angstrom", 20.0, float, 1.0, maximum=MAX_BOX_SIZE)
        if s.anchor_mode != "spread":  # spread has no random element
            say("\n  A random seed makes the anchor choice repeatable. Leave blank for none.")
            s.seed = ask_number("Random seed (0 or higher)", None, int, 0, required=False)
        if ask_yes_no("Show expert options (manual core pattern, debug output)?", False):
            expert_options(s)
    if s.size <= s.extent:
        warn(f"The search box ({s.size} A) is not larger than the template ({s.extent:.1f} A at its longest). "
             "Use a bigger box.")
    if s.size >= LARGE_BOX_SIZE:
        warn(f"A {s.size} A box is large: grid generation time and disk use grow with the cube of the size "
             f"(about {(s.size / 20.0) ** 3:.0f}x the default 20 A box).")
    many = len(s.kept) >= 10
    per_ligand = disk_per_ligand_mb(s.size, anchors_in_use(s))
    say(f"\n  Intermediate files (grid maps etc.) take about {per_ligand:.0f} MB per ligand at this box size.")
    say("  Deleting them keeps only the docked poses, the template and the receptor.")
    s.clean = ask_yes_no("Delete intermediate files when each ligand is done?", many)


def expert_options(s):
    say("\n  The core pattern must match the template and be valid SMILES. Leave it blank to find the core")
    say("  automatically.")
    available = []

    def validate_pattern(text):
        available[:] = manual_core_anchors(text, s.template_mol)
        return text

    s.query_smarts = ask("Manual core pattern (SMARTS)", validate=validate_pattern, required=False)
    if s.query_smarts:
        say(f"  Atom positions in that pattern that can be anchors: {sorted(available)}")
        text = ask("Anchor atom positions, space-separated (0-based), blank for automatic",
                   validate=lambda t: " ".join(map(str, parse_anchor_indices(t, set(available)))), required=False)
        s.anchor_indices = [int(t) for t in text.split()] if text else None
        warn("Manual pattern: the ligand checks above used automatic core finding and may not apply.")
    s.verbose = ask_yes_no("Print extra debug information?", False)


def autodock_gpu_text(s):
    """Which AutoDock-GPU the job will use, in words."""
    if s.adg_path:
        return f"custom build in {s.adg_path} ({atom_type_limit(s)} atom types)"
    return f"cluster module ({MAX_ATOM_TYPES} atom types)"


def cluster_settings_rows(s):
    """The settings that the 'change partition, GPUs, modules or pixi location' question is about, as they stand."""
    modules = ", ".join(m.strip() for m in s.modules.split(",") if m.strip()) or "none"
    return [("Partition", s.partition),
            ("Reservation", s.reservation or "none"),
            ("GPUs", s.gpus),
            ("Modules", modules),
            ("pixi folder", s.pixi_dir or "not found on PATH (you will be asked for it)"),
            ("AutoDock-GPU", autodock_gpu_text(s))]


def step_run_mode(s):
    heading(5, "How to run it")
    has_slurm = shutil.which("sbatch") is not None
    default_name = "dock_" + re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.splitext(os.path.basename(s.ligands))[0])
    s.jobname = default_name[:40]
    options = []
    if has_slurm:
        options.append(("submit", "Submit to the cluster queue (SLURM)"))
    options += [("run", "Run now in this terminal (only on a machine with a GPU and AutoDock-GPU loaded)"),
                ("save", "Only save the job script, do not run anything")]
    say("  Where should the docking run?")
    for n, (_, text) in enumerate(options, 1):
        say(f"    {n}. {text}")
    if not has_slurm:
        warn("SLURM (sbatch) was not found on this computer, so cluster submission is unavailable.")
    def pick(text):
        if text.isdigit() and 1 <= int(text) <= len(options):
            return options[int(text) - 1][0]
        raise ValueError(f"Please enter a number from 1 to {len(options)}.")

    s.mode = ask("Choose a number", "1", pick)

    say("\n  Results are written to an output folder; one numbered sub-folder per ligand.")
    protected = [s.receptor, s.template, s.ligands, REPO_DIR, os.getcwd(), os.path.expanduser("~"), "/"]
    if s.structure_path:
        protected.append(s.structure_path)
    default_out = os.path.join(os.path.dirname(s.ligands), "output")
    while True:
        s.resume, s.remaining = False, len(s.kept)
        out = ask("Output folder", default_out,
                  lambda t: os.path.abspath(os.path.expanduser(t.strip("'\""))))
        problem = validate_outdir(out, protected)
        if problem:
            fail(problem)
            continue
        if os.path.isdir(out) and os.listdir(out):
            warn(f"{out} already exists and is not empty.", record=False)
            earlier = previous_run(s, out)
            if earlier is not None:
                differences, done = earlier
                if differences:
                    say("       It holds an earlier run, but these differ from your current choices, so it "
                        f"cannot be resumed: {', '.join(differences)}.")
                else:
                    say(f"       It holds an earlier run with the same inputs and settings: {done} of "
                        f"{len(s.kept)} ligand(s) are already docked.")
                    if ask_yes_no("Resume it? Finished ligands are skipped, the rest are docked", True):
                        s.resume, s.force, s.outdir = True, False, out
                        s.remaining = len(s.kept) - done
                        break
            if not ask_yes_no("Delete its contents and write new results there?", False):
                default_out = out.rstrip("/") + "_new"
                continue
        s.force = os.path.isdir(out)
        s.outdir = out
        break

    needed = estimate_disk_mb(s.remaining, s.size, s.clean, anchors_in_use(s))
    free = shutil.disk_usage(nearest_existing_parent(s.outdir)).free / 1e6
    if needed > free:
        warn(f"This run may need about {needed / 1000:.1f} GB but only {free / 1000:.1f} GB is free where the "
             "output goes. Free up space, choose another folder, or delete intermediate files.")
    elif needed > 0.5 * free:
        warn(f"This run may need about {needed / 1000:.1f} GB, over half of the {free / 1000:.1f} GB free "
             "where the output goes.")

    if s.mode == "run":
        return
    s.time_limit = fmt_minutes(max(MIN_TIME_MINUTES, s.remaining * MINUTES_PER_LIGAND))
    say("\n  Cluster job settings (press Enter to accept the suggestion).")
    s.jobname = ask("Job name", s.jobname, lambda t: re.sub(r"[^A-Za-z0-9_.-]", "_", t))
    say(f"  Time limit: budgeting {MINUTES_PER_LIGAND} minutes per ligand (the example needs ~2).")
    s.time_limit = ask("Time limit (hours:minutes:seconds)", s.time_limit, time_validator)
    pixi = shutil.which("pixi")
    s.pixi_dir = os.path.dirname(pixi) if pixi else ""
    say("\n  The job will use these cluster settings:")
    table(cluster_settings_rows(s), indent="    ")
    if ask_yes_no("Change partition, GPUs, modules or pixi location?", False):
        s.partition = ask("Partition", s.partition)
        s.reservation = ask_reservation(s)
        s.gpus = ask_number("Number of GPUs", s.gpus, int, 1)
        s.modules = ask("Modules to load (comma-separated; GCC must come first)", s.modules)
        s.pixi_dir = ask("Folder containing the pixi program", s.pixi_dir or os.path.expanduser("~/.local/bin"))
        if s.adg_path and not ask_yes_no("Use the AutoDock-GPU with the higher atom-type limit that some ligands need?", True):
            s.adg_path = ""
    elif not s.pixi_dir:
        s.pixi_dir = ask("Folder containing the pixi program", os.path.expanduser("~/.local/bin"))
    ensure_partition_access(s)
    problem = build_problem(s)
    if problem:
        warn(problem, record=False)
        if ask_yes_no("Use the cluster's AutoDock-GPU instead (a ligand may then get fewer anchors)?", True):
            s.adg_path = ""


def output_note(s):
    if s.resume:
        return f" (resuming: {len(s.kept) - s.remaining} finished ligand(s) are skipped)"
    if s.force and os.listdir(s.outdir):
        return " (existing contents will be deleted)"
    return ""


def anchors_text(s):
    if s.anchor_mode == "spread":
        return f"up to {s.anchor_num} per ligand, spread across the shared core"
    return f"{s.anchor_num} ({s.anchor_mode})"


def step_review(s):
    heading(6, "Review")
    problem = build_problem(s)
    if problem:
        warn(problem)
    names = ", ".join(r.name for r in s.kept[:3]) + (" ..." if len(s.kept) > 3 else "")
    if s.prepared is not None:
        receptor_path, template_path, _ = prepared_paths(s)
        rows = [("Structure", s.structure_path),
                ("Receptor", f"{receptor_path} (chains {', '.join(s.prepared.kept_chains)}, cleaned)"),
                ("Template", f"{template_path} ({s.prepared.template.label})")]
    elif s.protonated is not None:
        rows = [("Receptor", f"{prepared_paths(s)[0]} (protonated from {s.original_receptor})"),
                ("Template", s.template)]
    else:
        rows = [("Receptor", s.receptor), ("Template", s.template)]
    rows += [("Ligands", f"{len(s.kept)} of {s.all_count} in {s.ligands} ({names})"),
            ("Anchors", anchors_text(s)), ("Search box", f"{s.size} A"),
            ("Output folder", s.outdir + output_note(s))]
    if s.seed is not None:
        rows.append(("Random seed", s.seed))
    rows.append(("AutoDock-GPU", autodock_gpu_text(s)))
    if s.protonated is not None:
        rows.append(("Protonation", f"PDB2PQR + PROPKA at pH {s.protonated.ph:g} ({s.protonated.hydrogens} hydrogens)"))
    elif s.add_hydrogens:
        rows.append(("Hydrogens", "added to the receptor by ADFR (quick rule-based guess)"))
    rows.append(("Delete intermediates", "yes" if s.clean else "no"))
    if s.mode != "run":
        rows += [("Partition / GPUs", f"{s.partition} / {s.gpus}" + (f" (reservation {s.reservation})" if s.reservation else "")),
                 ("Time limit", s.time_limit), ("Job name", s.jobname)]
    action = {"submit": "Submit the job to the queue", "run": "Run the docking now",
              "save": "Save the job script"}[s.mode]
    rows.append(("Action", action))
    table(rows, indent="  ")
    if len(s.kept) != s.all_count:
        mapping = ", ".join(f"{n} = {r.name[:25]}" for n, r in enumerate(s.kept[:8]))
        say("\n  Some ligands are left out, so result folders are numbered within the reduced list:")
        say(f"    {mapping}{' ...' if len(s.kept) > 8 else ''}")
        say("  The reduced list is saved as <job name>_ligands.sdf next to the job script.")
    if WARNINGS:
        say(f"\n  {_c('33', 'Warnings to consider:')}")
        for w in dict.fromkeys(WARNINGS):
            say(f"   - {w}")
    say("")
    return ask_yes_no(f"{action}?", True)


# --------------------------------------------------------------------------
# Running, submitting, watching
# --------------------------------------------------------------------------


def write_script(s):
    finalize_structure_files(s)
    prepare_ligand_file(s)
    os.makedirs(jobdir(s), exist_ok=True)
    s.script_path = os.path.join(jobdir(s), f"submit_{s.jobname}.sh")
    with open(s.script_path, "w") as fh:
        fh.write(build_script(s))
    os.chmod(s.script_path, 0o755)
    ok(f"Job script saved: {s.script_path}")


def latest_job_id(jobname):
    out = subprocess.run(["squeue", "-h", "-u", getpass.getuser(), "-n", jobname, "-o", "%i|%V"],
                         capture_output=True, text=True).stdout
    jobs = [line.split("|") for line in out.splitlines() if "|" in line]
    return max(jobs, key=lambda j: j[1])[0] if jobs else None


def submit(s):
    say("")
    result = subprocess.run(["sbatch", s.script_path])  # inherits the terminal: sbatch may ask questions
    if result.returncode != 0:
        fail("sbatch reported an error, so the job was not submitted.")
        why = check_submission(s.partition, s.reservation, s.gpus, s.time_limit)
        if why.ok is False:
            say(f"  Reason: {why.reason}. Run the assistant again and choose a partition you can use.")
        return None
    job_id = latest_job_id(s.jobname)
    if job_id is None:
        warn("Could not find the job in the queue. If you declined at the submit prompt, it was not submitted.")
    return job_id


def job_state(job_id):
    """Ask SLURM about a job. Returns ("state", [state, run time]), ("gone", None) or ("error", message).

    "gone" means SLURM answered and the job is no longer in the queue. A failure to
    reach SLURM is an "error", never "gone", so a hiccup is not mistaken for a finished job.
    """
    try:
        result = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%T|%M"],
                                capture_output=True, text=True, check=False)
    except OSError as exc:
        return "error", str(exc)
    if result.returncode == 0:
        out = result.stdout.strip()
        return ("state", out.split("|")) if out else ("gone", None)
    if "invalid job id" in result.stderr.lower():
        return "gone", None
    return "error", result.stderr.strip() or f"squeue exited with code {result.returncode}"


def job_final_state(job_id):
    """(state, exit code) from the accounting database, or None if it is unavailable."""
    try:
        result = subprocess.run(["sacct", "-n", "-X", "-P", "-j", job_id, "-o", "State,ExitCode"],
                                capture_output=True, text=True, check=False)
    except OSError:
        return None
    lines = result.stdout.strip().splitlines() if result.returncode == 0 else []
    return tuple(lines[0].split("|")[:2]) if lines and "|" in lines[0] else None


def watch(job_id, s):
    say(f"\n  Watching job {job_id}. Press Ctrl+C to stop watching (the job keeps running).\n")
    last_state, last_print, errors = None, 0.0, 0
    try:
        while True:
            kind, data = job_state(job_id)
            if kind == "error":
                errors += 1
                if errors >= MAX_WATCH_ERRORS:
                    warn(f"Could not reach SLURM to check the job ({data}). Stopped watching; "
                         "your job is not affected.", record=False)
                    show_follow_up(job_id, s)
                    return False
                time.sleep(POLL_SECONDS)
                continue
            errors = 0
            if kind == "state":
                state = data
                now = time.time()
                if state[0] != last_state or now - last_print > 60:
                    label = {"PENDING": "waiting in the queue", "RUNNING": "running"}.get(state[0], state[0].lower())
                    say(f"  {time.strftime('%H:%M:%S')}  {label} (run time {state[1]})")
                    last_state, last_print = state[0], now
                time.sleep(POLL_SECONDS)
                continue
            final = job_final_state(job_id)  # "gone" from the queue; make sure it is really over
            if final and final[0].split()[0] in ACTIVE_STATES:
                time.sleep(POLL_SECONDS)
                continue
            break
    except KeyboardInterrupt:
        say("\n  Stopped watching. The job is still in the queue or running.")
        show_follow_up(job_id, s)
        return False
    if final:
        say(f"\n  Job finished: {final[0]} (exit code {final[1]})")
    else:
        say("\n  Job finished.")
    return True


def show_follow_up(job_id, s):
    log = f"{os.path.join(jobdir(s), s.jobname)}_{job_id}.log"
    say("\n  Useful commands:")
    say(f"    squeue -u {getpass.getuser()}            show your queued and running jobs")
    say(f"    tail -f {log}   follow the job's log (Ctrl+C to stop)")
    say(f"    scancel {job_id}                         cancel the job")
    if mail_address():
        say(f"\n  An email goes to {mail_address()} when the job ends or fails.")


def result_rows(s):
    """One row per ligand: from summary.tsv written by the job, falling back to the files on disk."""
    import anchored_docking as ad
    import results

    summary = {int(r["index"]): r for r in ad.read_summary(s.outdir) if r.get("index", "").isdigit()}
    rows = []
    for n, ligand in enumerate(s.kept):
        info = summary.get(n)
        pose_file = os.path.join(s.outdir, str(n), "ligand_docked.sdf")
        has_pose = os.path.isfile(pose_file)
        energy = results.top_pose_scores(pose_file).get("free_energy", "") if has_pose else ""
        if info is None:
            status = "docked (not checked)" if has_pose else "not run"
            rows.append({"index": n, "name": ligand.name, "status": status, "energy": energy, "anchor": "",
                         "rmsd": "", "note": ""})
            continue
        status = {"docked": "docked", "flagged": "CHECK POSE", "failed": "FAILED"}.get(info["status"], info["status"])
        rows.append({"index": n, "name": ligand.name, "status": status, "energy": energy,
                     "anchor": info.get("anchor_max_dev", ""),
                     "rmsd": info.get("core_rmsd", ""), "note": info.get("note", "")})
    return rows


def show_results(s, job_id=None):
    rows = result_rows(s)
    say("\n  Results:")
    table([(r["index"], r["name"][:30], r["status"], r["energy"], r["anchor"], r["rmsd"]) for r in rows],
          headers=("Folder", "Ligand", "Result", "Energy (kcal/mol)", "Anchor dev (A)", "Core RMSD (A)"))
    docked = [r for r in rows if r["status"] in ("docked", "CHECK POSE", "docked (not checked)")]
    flagged = [r for r in rows if r["status"] == "CHECK POSE"]
    failed = [r for r in rows if r["status"] == "FAILED"]
    done = len(docked)
    say(f"\n  {done} of {len(s.kept)} ligand(s) have a docked pose in {s.outdir}")
    if flagged:
        say(f"\n  {_c('33', str(len(flagged)) + ' pose(s) need a closer look.')} Docking only pulls the anchors "
            "towards the template, so a pose can still end up elsewhere:")
        for r in flagged[:10]:
            say(f"    folder {r['index']} ({r['name'][:30]}): {r['note']}")
    if failed:
        say(f"\n  {_c('31', str(len(failed)) + ' ligand(s) failed.')} The other ligands were still docked:")
        for r in failed[:10]:
            say(f"    folder {r['index']} ({r['name'][:30]}): {r['note']}")
        say(f"  Details are in error.txt inside each failed folder. After fixing the cause, run `{RUN_COMMAND}` "
            "again with the same folder and answer yes to resume, to redo only the failed ligands.")
    if done:
        say(f"\n  Files in {s.outdir}:")
        say("    results.csv                 one row per ligand: status, energies, pose check")
        say("    top_poses.sdf               the best pose of every docked ligand, with scores (all hydrogens)")
        say("    all_poses.sdf               every cluster-lead pose of every ligand")
        say("    view_all.pml                opens the protein, template and all ligands in one PyMOL session")
        say("    <folder>/ligand_docked.sdf  the poses of one ligand;  <folder>/view_results.pml  to view just that one")
        say("  The protein and template used for viewing have every hydrogen back (AutoDock merges the "
            "non-polar ones away).")
        say("  Energies include the anchoring bias, so use them to rank within this run, not against unbiased docking.")
    if done < len(s.kept) and job_id:
        log = f"{os.path.join(jobdir(s), s.jobname)}_{job_id}.log"
        if os.path.isfile(log):
            say(f"\n  Last lines of {log}:")
            with open(log, errors="replace") as fh:
                for line in fh.readlines()[-15:]:
                    say("    " + line.rstrip())


def run_here(s):
    sys.path.insert(0, REPO_DIR)
    import anchored_docking as ad

    try:
        ad.find_autodock_gpu()
    except FileNotFoundError as exc:
        fail(str(exc))
        return None
    if shutil.which("nvidia-smi") is None or subprocess.run(["nvidia-smi", "-L"], capture_output=True).returncode != 0:
        warn("No GPU is visible on this machine; AutoDock-GPU will fail here.")
        if not ask_yes_no("Run anyway?", False):
            return None
    finalize_structure_files(s)
    prepare_ligand_file(s)
    os.makedirs(jobdir(s), exist_ok=True)
    say("\n  Running... (Ctrl+C to stop)\n")
    try:
        result = subprocess.run([sys.executable, os.path.join(REPO_DIR, "anchored_docking.py")] + build_args(s),
                                cwd=REPO_DIR)
    except KeyboardInterrupt:
        fail("Stopped.")
        return None
    return result.returncode == 0


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def check_environment():
    try:
        import gemmi  # noqa: F401
        import meeko  # noqa: F401
        import rdkit  # noqa: F401
        import scipy  # noqa: F401
    except ImportError:
        print("This tool must be started through pixi so that RDKit and Meeko are available:\n\n"
              f"    {RUN_COMMAND}\n")
        sys.exit(1)
    if not adfr_ready():
        if UNDER_DDSUITE:
            print("ADFRsuite, which prepares the receptor, is not installed in this copy of the tool. "
                  "Ask the maintainer of the suite to reinstall it.")
        else:
            print("ADFRsuite, which prepares the receptor, is not installed yet. Install everything once with:\n\n"
                  "    pixi install --all\n\n"
                  f"then start this tool again with `{RUN_COMMAND}`.")
        sys.exit(1)


def main():
    argparse.ArgumentParser(
        description=f"Interactive anchored docking. Run it with `{RUN_COMMAND}` and answer the questions."
    ).parse_args()
    check_environment()

    s = Settings()
    say(bold("\nAnchored docking"))
    say("Dock a series of ligands into a protein, anchored to a template ligand.")
    say("Answer each question, or press Enter to accept the suggestion in [brackets]. Ctrl+C quits.")
    s_submitted = False
    try:
        step_files(s)
        step_settings(s)
        step_run_mode(s)
        if not step_review(s):
            say("\nCancelled. Nothing was run.")
            return
        s_submitted = True
        if s.mode == "run":
            outcome = run_here(s)
            if outcome is None:
                fail("The docking did not start or was stopped. See the messages above.")
                sys.exit(1)
            show_results(s)  # also after failures: the other ligands may have been docked
            if not outcome:
                fail("Some ligands failed. See the results above.")
                sys.exit(1)
            return
        write_script(s)
        if s.mode == "save":
            say(f"\n  Submit it later with:  sbatch {s.script_path}")
            return
        job_id = submit(s)
        if job_id is None:
            return
        ok(f"Submitted job {job_id}.")
        show_follow_up(job_id, s)
        say("")
        if ask_yes_no("Watch the job until it finishes?", True):
            if watch(job_id, s):
                show_results(s, job_id)
    except (KeyboardInterrupt, EOFError):
        say("\n\nStopped." + ("" if s_submitted else " Nothing was run."))
        sys.exit(130)


if __name__ == "__main__":
    main()
