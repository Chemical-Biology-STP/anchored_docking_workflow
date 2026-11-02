"""Check docked poses against the template and the anchors they were meant to satisfy.

Docking only pulls the anchor atoms towards their targets through a soft potential
(scripts/addbias.py), so a pose can come out with the core in the wrong place and
still be written as a normal result. This module measures how far a pose is from
what the anchoring asked for, so that such poses can be flagged.
"""

import numpy as np
from rdkit import Chem

# scripts/addbias.py adds no penalty within 1.2 A of an anchor target and a steep one (100 kcal/mol
# per A) beyond it, but the 0.375 A grid interpolation blurs that edge. A pose that satisfies the
# anchoring therefore ends up at up to about 1.3 A: two real AutoDock-GPU poses of the example
# ligand measured 1.32 A and 1.24 A. 2.0 A leaves a margin above that, so only poses that clearly
# escaped the anchoring are flagged. This is a judgement based on those two runs, not a calibration.
ANCHOR_TOLERANCE = 2.0
# Heuristic: RMSD over all matched core atoms (the ligand core is not identical to the template's).
CORE_RMSD_TOLERANCE = 2.5
MAX_MATCHES = 5000


def read_poses(path):
    """All poses in a docked SDF file (hydrogens kept). Unreadable entries are dropped."""
    try:
        return [m for m in Chem.SDMolSupplier(path, removeHs=False) if m is not None]
    except OSError:  # RDKit refuses empty or unreadable files
        return []


def core_matches(pose, pattern_mol, smarts):
    """Ways to map the ligand core (atoms of pattern_mol) onto a pose, as dicts pattern atom -> pose atom.

    Uses the same pattern matching as anchored_docking.find_anchor_points, so the mapping
    agrees with the one used to pick the anchors. Symmetric cores give several mappings.
    """
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return []
    pattern_order = pattern_mol.GetSubstructMatch(query)
    if not pattern_order:
        return []
    matches = pose.GetSubstructMatches(query, uniquify=False, maxMatches=MAX_MATCHES)
    return [{p: m[i] for i, p in enumerate(pattern_order)} for m in matches]


def check_pose(pose, template_mol, pattern_mol, smarts, ptn2template, anchor_indices, anchor_positions,
               anchor_tolerance=ANCHOR_TOLERANCE, core_tolerance=CORE_RMSD_TOLERANCE):
    """Measure one pose.

    anchor_indices are pattern-atom indices and anchor_positions the target coordinates
    (taken from the template). Returns a dict with the numbers and a verdict.
    When the core is symmetric, the mapping that fits the anchors best is used, so the
    result is the most favourable reading of the pose.
    """
    result = {"anchor_deviations": [], "anchor_max": None, "core_rmsd": None, "core_atoms": 0,
              "ok": None, "problems": []}
    mappings = core_matches(pose, pattern_mol, smarts)
    if not mappings:
        result["problems"].append("could not map the ligand core onto the docked pose, so it was not checked")
        return result

    pose_xyz = pose.GetConformer().GetPositions()
    template_xyz = template_mol.GetConformer().GetPositions()
    targets = np.array(anchor_positions, dtype=float).reshape(-1, 3)

    best = None
    for mapping in mappings:
        deviations = [float(np.linalg.norm(pose_xyz[mapping[a]] - t))
                      for a, t in zip(anchor_indices, targets) if a in mapping]
        pairs = [(pose_xyz[mapping[p]], template_xyz[t]) for p, t in ptn2template.items() if p in mapping]
        rmsd = (float(np.sqrt(np.mean([np.sum((a - b) ** 2) for a, b in pairs]))) if pairs else None)
        key = (max(deviations) if deviations else 0.0, rmsd if rmsd is not None else 0.0)
        if best is None or key < best[0]:
            best = (key, deviations, rmsd, len(pairs))
    _, deviations, rmsd, n_atoms = best

    result["anchor_deviations"] = [round(x, 2) for x in deviations]
    result["anchor_max"] = round(max(deviations), 2) if deviations else None
    result["core_rmsd"] = round(rmsd, 2) if rmsd is not None else None
    result["core_atoms"] = n_atoms
    if result["anchor_max"] is not None and result["anchor_max"] > anchor_tolerance:
        result["problems"].append(f"an anchor atom is {result['anchor_max']} A from its target "
                                  f"(expected at most {anchor_tolerance} A)")
    if rmsd is not None and rmsd > core_tolerance:
        result["problems"].append(f"the core is {result['core_rmsd']} A RMSD from the template "
                                  f"(expected at most {core_tolerance} A)")
    result["ok"] = not result["problems"]
    return result


def check_docked_file(path, template_mol, pattern_mol, smarts, ptn2template, anchor_indices, anchor_positions):
    """Check every pose in a docked SDF. The first pose is the top-ranked one."""
    poses = read_poses(path)
    if not poses:
        return {"ok": None, "problems": ["the docked file has no readable pose"], "poses": 0,
                "anchor_deviations": [], "anchor_max": None, "core_rmsd": None, "core_atoms": 0,
                "poses_ok": 0}
    checks = [check_pose(p, template_mol, pattern_mol, smarts, ptn2template, anchor_indices, anchor_positions)
              for p in poses]
    first = dict(checks[0])
    first["poses"] = len(checks)
    first["poses_ok"] = sum(1 for c in checks if c["ok"])
    return first
