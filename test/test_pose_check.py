"""Tests for pose_check.py: does a docked pose satisfy the anchoring it was meant to?"""

import os

import numpy as np
import pytest
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import anchored_docking as ad
import pose_check as pc

RDLogger.DisableLog("rdApp.*")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(REPO, "test", "crystal_ligand.sdf")
LIGAND = os.path.join(REPO, "test", "CAT-13f.sdf")


@pytest.fixture(scope="module")
def template():
    return Chem.SDMolSupplier(TEMPLATE)[0]


@pytest.fixture(scope="module")
def search(template):
    ligand = Chem.SDMolSupplier(LIGAND)[0]
    return ad.find_anchor_points(template, ligand, 3, 100, 1)


@pytest.fixture(scope="module")
def good_pose():
    """The example ligand in the crystal frame: a stand-in for a pose that satisfies the anchoring."""
    return Chem.SDMolSupplier(LIGAND, removeHs=False)[0]


def anchors(search):
    return [idx for _, idx in search.anchor_patterns]


def run_check(pose, template, search, **kwargs):
    return pc.check_pose(pose, template, search.pattern_mol, search.pattern_smarts, search.ptn2template,
                         anchors(search), search.anchor_positions, **kwargs)


def transformed(mol, matrix=None, shift=(0.0, 0.0, 0.0)):
    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    xyz = conf.GetPositions()
    centre = xyz.mean(axis=0)
    for i, p in enumerate(xyz):
        q = (np.asarray(matrix) @ (p - centre) + centre if matrix is not None else p) + np.asarray(shift)
        conf.SetAtomPosition(i, tuple(q))
    return mol


def test_anchor_search_exposes_the_core_mapping(search, template):
    assert len(search.anchor_patterns) == 3 and len(search.anchor_positions) == 3
    assert search.ptn2template and all(0 <= t < template.GetNumAtoms() for t in search.ptn2template.values())
    assert search.pattern_mol.GetNumHeavyAtoms() == len(search.ptn2template)
    assert Chem.MolFromSmarts(search.pattern_smarts) is not None


def test_pose_on_the_anchors_passes(good_pose, template, search):
    result = run_check(good_pose, template, search)
    assert result["ok"] is True and result["problems"] == []
    assert result["anchor_max"] < pc.ANCHOR_TOLERANCE
    assert result["core_rmsd"] < pc.CORE_RMSD_TOLERANCE
    assert result["core_atoms"] == search.pattern_mol.GetNumHeavyAtoms()
    assert len(result["anchor_deviations"]) == 3


def test_shifted_pose_is_flagged(good_pose, template, search):
    result = run_check(transformed(good_pose, shift=(6.0, 0.0, 0.0)), template, search)
    assert result["ok"] is False
    assert result["anchor_max"] == pytest.approx(6.0, abs=0.7)
    assert any("anchor atom" in p for p in result["problems"])
    assert any("RMSD" in p for p in result["problems"])


def test_flipped_pose_is_flagged(good_pose, template, search):
    flip = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    result = run_check(transformed(good_pose, matrix=flip), template, search)
    assert result["ok"] is False and result["core_rmsd"] > pc.CORE_RMSD_TOLERANCE


def test_tolerances_can_be_tightened(good_pose, template, search):
    result = run_check(good_pose, template, search, anchor_tolerance=0.1, core_tolerance=0.1)
    assert result["ok"] is False and len(result["problems"]) == 2


def test_unrelated_pose_is_reported_as_not_checked(template, search):
    ethanol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(ethanol, AllChem.ETKDGv3())
    result = run_check(ethanol, template, search)
    assert result["ok"] is None and "could not map" in result["problems"][0]


def test_check_works_with_explicit_hydrogens(good_pose, template, search):
    with_h = Chem.AddHs(good_pose, addCoords=True)
    assert with_h.GetNumAtoms() > good_pose.GetNumAtoms()
    assert run_check(with_h, template, search)["ok"] is True


def test_docked_file_reports_first_pose_and_counts_good_poses(good_pose, template, search, tmp_path):
    path = tmp_path / "docked.sdf"
    writer = Chem.SDWriter(str(path))
    writer.write(good_pose)
    writer.write(transformed(good_pose, shift=(8.0, 0.0, 0.0)))
    writer.close()
    result = pc.check_docked_file(str(path), template, search.pattern_mol, search.pattern_smarts,
                                  search.ptn2template, anchors(search), search.anchor_positions)
    assert result["poses"] == 2 and result["poses_ok"] == 1 and result["ok"] is True


def test_empty_docked_file_is_not_a_pass(template, search, tmp_path):
    path = tmp_path / "empty.sdf"
    path.write_text("")
    result = pc.check_docked_file(str(path), template, search.pattern_mol, search.pattern_smarts,
                                  search.ptn2template, anchors(search), search.anchor_positions)
    assert result["ok"] is None and result["poses"] == 0


def test_docking_script_saves_the_check_and_summary_row(good_pose, template, search, tmp_path):
    docked = tmp_path / "ligand_docked.sdf"
    writer = Chem.SDWriter(str(docked))
    writer.write(transformed(good_pose, shift=(7.0, 0.0, 0.0)))
    writer.close()
    checked = ad.check_docked_pose(str(tmp_path), template, search)
    assert checked["ok"] is False and (tmp_path / "pose_check.json").is_file()
    row = ad.row_from_check(0, "ligand", checked)
    assert row["status"] == "flagged" and row["anchor_max_dev"] > pc.ANCHOR_TOLERANCE
    assert "anchor atom" in row["note"]
