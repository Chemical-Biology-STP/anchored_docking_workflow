"""Tests for results.py (combined outputs and display files with every hydrogen) and the template
hydrogens taken from the PDB's component dictionary (structure_prep.py).

Nothing here uses the network: dictionary entries are built in the tests.
"""

import csv
import json
import os
import subprocess

import gemmi
import numpy as np
import pytest
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import anchored_docking as ad
import results
import structure_prep as sp

RDLogger.DisableLog("rdApp.*")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(REPO, "test", "Bace1_protein.pdb")
TEMPLATE = os.path.join(REPO, "test", "crystal_ligand.sdf")
LIGAND = os.path.join(REPO, "test", "CAT-13f.sdf")


# --------------------------------------------------------------------------
# Protein: PDBQT -> PDB and the receptor with every hydrogen
# --------------------------------------------------------------------------

PDBQT = (
    "ATOM      1  N   ALA A   1      10.000  11.000  12.000  1.00  0.00    -0.350 N \n"
    "ATOM      2  CA  ALA A   1      11.000  11.500  12.500  1.00  0.00     0.180 C \n"
    "ATOM      3  CG  PHE A   2      12.000  12.000  13.000  1.00  0.00    -0.020 A \n"
    "ATOM      4  OG  SER A   3      13.000  13.000  14.000  1.00  0.00    -0.390 OA\n"
    "ATOM      5  HG  SER A   3      13.500  13.500  14.500  1.00  0.00     0.210 HD\n"
    "ATOM      6  HA  ALA A   1      11.200  11.900  12.100  1.00  0.00     0.000 H \n"
    "ATOM      7  SD  MET A   4      14.000  14.000  15.000  1.00  0.00    -0.090 SA\n"
    "ATOM      8  ND1 HIS A   5      15.000  15.000  16.000  1.00  0.00    -0.250 NA\n"
    "HETATM    9 ZN    ZN A 301      16.000  16.000  17.000  1.00  0.00     2.000 Zn\n"
    "HETATM   10 CL    CL A 302      17.000  17.000  18.000  1.00  0.00    -1.000 Cl\n"
)


def test_pdbqt_types_become_elements():
    pdb = results.pdbqt_to_pdb(PDBQT)
    elements = [l[76:78].strip() for l in pdb.splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert elements == ["N", "C", "C", "O", "H", "H", "S", "N", "ZN", "CL"]
    first = pdb.splitlines()[0]
    assert first[12:16].strip() == "N" and first[17:20] == "ALA" and float(first[30:38]) == 10.0
    assert "-0.350" not in pdb  # the AutoDock charge and type columns are gone
    assert pdb.rstrip().endswith("END")
    assert pdb.splitlines()[8].startswith("HETATM")


def test_converted_receptor_is_readable_and_keeps_hydrogens():
    pdb = results.pdbqt_to_pdb(PDBQT)
    mol = Chem.MolFromPDBBlock(pdb, sanitize=False, removeHs=False)
    assert mol.GetNumAtoms() == 10
    assert sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) == 2


def test_hydrogens_are_counted_by_element_or_by_name(tmp_path):
    path = tmp_path / "r.pdb"
    path.write_text(
        "ATOM      1  N   ALA A   1      10.000  11.000  12.000  1.00  0.00           N\n"
        "ATOM      2  H   ALA A   1      10.500  11.500  12.500  1.00  0.00           H\n"
        "ATOM      3  HA  ALA A   1      10.500  11.500  12.500  1.00  0.00\n"
        "ATOM      4  CA  ALA A   1      11.000  11.000  12.000  1.00  0.00\n")
    assert results.count_hydrogens(str(path)) == 2


def test_a_receptor_that_has_hydrogens_is_used_as_it_is(tmp_path, monkeypatch):
    monkeypatch.setattr(results.subprocess, "run", lambda *a, **k: pytest.fail("ADFR must not run"))
    destination = tmp_path / "view.pdb"
    n, how = results.make_view_receptor(RECEPTOR, "/adfr", True, str(destination))
    assert n == 2975 and "as supplied" in how
    assert destination.read_text() == open(RECEPTOR).read()


def stripped_receptor(tmp_path):
    path = tmp_path / "noH.pdb"
    with open(RECEPTOR) as src, open(path, "w") as dst:
        dst.writelines(l for l in src if l.startswith("ATOM") and l[76:78].strip() != "H")
    return str(path)


def test_a_receptor_without_hydrogens_stays_without_when_none_were_added(tmp_path, monkeypatch):
    monkeypatch.setattr(results.subprocess, "run", lambda *a, **k: pytest.fail("ADFR must not run"))
    n, how = results.make_view_receptor(stripped_receptor(tmp_path), "/adfr", False, str(tmp_path / "v.pdb"))
    assert n == 0 and "none were added" in how


def test_added_hydrogens_come_from_a_run_that_does_not_merge_them(tmp_path, monkeypatch):
    calls = []

    def fake(command, cwd, **kwargs):
        calls.append(command)
        with open(os.path.join(cwd, "full.pdbqt"), "w") as fh:
            fh.write(PDBQT)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(results.subprocess, "run", fake)
    destination = tmp_path / "v.pdb"
    n, how = results.make_view_receptor(stripped_receptor(tmp_path), "/adfr", True, str(destination))
    command = calls[0]
    assert command[command.index("-A") + 1] == "checkhydrogens"
    assert command[command.index("-U") + 1] == "deleteAltB"  # not the default, which merges non-polar hydrogens
    assert n == 2 and "non-polar" in how
    assert results.count_hydrogens(str(destination)) == 2


@pytest.mark.skipif(not os.access(os.path.join(ad.DEFAULT_ADFR_PATH, "bin", "prepare_receptor"), os.X_OK),
                    reason="ADFRsuite is not installed (pixi install --all)")
def test_real_adfr_adds_all_hydrogens_and_keeps_the_docking_ones_in_place(tmp_path):
    source = stripped_receptor(tmp_path)
    destination = tmp_path / "view.pdb"
    n, _ = results.make_view_receptor(source, ad.DEFAULT_ADFR_PATH, True, str(destination))
    # what docking uses: the same receptor prepared the normal way, non-polar hydrogens merged away
    docking = tmp_path / "docking"
    docking.mkdir()
    ad.prepare_receptor(str(docking), source, ad.DEFAULT_ADFR_PATH, add_hydrogens=True)
    polar = {}
    for line in open(docking / "receptor.pdbqt"):
        if line.startswith("ATOM") and line.split()[-1] == "HD":
            polar[(line[12:16].strip(), line[17:20], line[21], line[22:26])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    shown = {}
    for line in open(destination):
        if line.startswith("ATOM") and line[76:78].strip() == "H":
            shown[(line[12:16].strip(), line[17:20], line[21], line[22:26])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    assert len(polar) > 400 and n > 2 * len(polar)            # many more hydrogens than the docking input has
    assert all(k in shown for k in polar)                      # every docking hydrogen is shown ...
    assert max(np.linalg.norm(polar[k] - shown[k]) for k in polar) < 0.001  # ... where docking had it


# --------------------------------------------------------------------------
# Template: hydrogens for display
# --------------------------------------------------------------------------


def test_template_without_hydrogens_gets_them_from_its_bond_orders(tmp_path):
    destination = tmp_path / "t.sdf"
    n, how = results.make_view_template(TEMPLATE, str(destination))
    mol = Chem.SDMolSupplier(str(destination), removeHs=False)[0]
    assert n == 19 and mol.GetNumAtoms() == 45 and "RDKit" in how


def test_template_that_already_has_hydrogens_is_not_changed(tmp_path):
    full = Chem.AddHs(Chem.SDMolSupplier(TEMPLATE)[0], addCoords=True)
    source = tmp_path / "full.sdf"
    with Chem.SDWriter(str(source)) as w:
        w.write(full)
    n, _ = results.make_view_template(str(source), str(tmp_path / "t.sdf"))
    assert n == 19 and Chem.SDMolSupplier(str(tmp_path / "t.sdf"), removeHs=False)[0].GetNumAtoms() == 45


def test_a_prepared_display_template_is_used_instead_of_guessing(tmp_path):
    prepared = tmp_path / "view.sdf"
    with Chem.SDWriter(str(prepared)) as w:
        w.write(Chem.AddHs(Chem.SDMolSupplier(TEMPLATE)[0], addCoords=True))
    n, how = results.make_view_template("unused.sdf", str(tmp_path / "t.sdf"), view_file=str(prepared))
    assert n == 19 and "prepared" in how


def test_a_display_template_without_hydrogens_is_reported_not_guessed(tmp_path):
    heavy = tmp_path / "heavy.sdf"
    heavy.write_text(Chem.MolToMolBlock(Chem.SDMolSupplier(TEMPLATE)[0]) + "$$$$\n")
    n, how = results.make_view_template("unused.sdf", str(tmp_path / "t.sdf"), view_file=str(heavy))
    assert n == 0 and "none could be derived" in how


# --------------------------------------------------------------------------
# Component dictionary
# --------------------------------------------------------------------------


def ccd_for(mol, drop_atom=None, rename=None, extra_h=0, comp="LIG"):
    """A component dictionary entry (gemmi block) describing `mol`, as the PDB would publish it."""
    full = Chem.AddHs(Chem.Mol(mol))
    kek = Chem.Mol(full)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    names = {}
    counts = {}
    for atom in kek.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol == "H":
            continue
        names[atom.GetIdx()] = f"{symbol}{atom.GetIdx() + 1}"
    h_i = 0
    for atom in kek.GetAtoms():
        if atom.GetSymbol() == "H":
            h_i += 1
            names[atom.GetIdx()] = f"H{h_i}"
    lines = ["data_" + comp, "loop_", "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
             "_chem_comp_atom.type_symbol", "_chem_comp_atom.charge"]
    for idx, name in names.items():
        a = kek.GetAtomWithIdx(idx)
        if name == drop_atom:
            continue
        lines.append(f"{comp} {rename.get(name, name) if rename else name} {a.GetSymbol()} {a.GetFormalCharge()}")
    for k in range(extra_h):
        lines.append(f"{comp} HX{k} H 0")
    lines += ["loop_", "_chem_comp_bond.comp_id", "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
              "_chem_comp_bond.value_order"]
    order = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB", Chem.BondType.TRIPLE: "TRIP"}
    for b in kek.GetBonds():
        lines.append(f"{comp} {names[b.GetBeginAtomIdx()]} {names[b.GetEndAtomIdx()]} {order[b.GetBondType()]}")
    return gemmi.cif.read_string("\n".join(lines) + "\n")[0]


def residue_mol(mol):
    """The template as it comes out of a structure file: heavy atoms, PDB names, single bonds."""
    pdb = Chem.MolToPDBBlock(mol)
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            idx = int(line[6:11]) - 1
            name = f"{mol.GetAtomWithIdx(idx).GetSymbol()}{idx + 1}"
            line = line[:12] + f" {name:<3}" + line[16:]
            out.append(line)
    return Chem.MolFromPDBBlock("\n".join(out) + "\n", sanitize=False, removeHs=False, proximityBonding=True)


@pytest.fixture(scope="module")
def real_template():
    return Chem.SDMolSupplier(TEMPLATE)[0]


def test_dictionary_chemistry_gives_the_right_hydrogens(real_template):
    result, why = sp.template_with_hydrogens(residue_mol(real_template), ccd_for(real_template))
    assert why == "" and result is not None
    assert sum(a.GetAtomicNum() == 1 for a in result.GetAtoms()) == 19
    assert result.GetNumHeavyAtoms() == real_template.GetNumHeavyAtoms()
    assert Chem.MolToSmiles(Chem.RemoveHs(result)) == Chem.MolToSmiles(real_template)  # real bond orders and charges
    assert any(a.GetIsAromatic() for a in result.GetAtoms())


def test_heavy_atom_coordinates_are_untouched_by_adding_hydrogens(real_template):
    mol = residue_mol(real_template)
    before = np.array(mol.GetConformer().GetPositions())
    result, _ = sp.template_with_hydrogens(mol, ccd_for(real_template))
    assert np.abs(result.GetConformer().GetPositions()[:mol.GetNumAtoms()] - before).max() < 1e-6


def test_hydrogens_sit_at_sensible_distances(real_template):
    result, _ = sp.template_with_hydrogens(residue_mol(real_template), ccd_for(real_template))
    xyz = result.GetConformer().GetPositions()
    lengths = [np.linalg.norm(xyz[b.GetBeginAtomIdx()] - xyz[b.GetEndAtomIdx()]) for b in result.GetBonds()
               if 1 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())]
    assert 0.95 < min(lengths) and max(lengths) < 1.15


def test_missing_ligand_atoms_are_not_papered_over(real_template):
    block = ccd_for(real_template)
    mol = residue_mol(real_template)
    editable = Chem.RWMol(mol)
    editable.RemoveAtom(0)  # the crystal lost an atom
    result, why = sp.template_with_hydrogens(editable.GetMol(), block)
    assert result is None and "does not match" in why and "missing" in why


def test_a_different_ligand_is_rejected(real_template):
    other = ccd_for(Chem.MolFromSmiles("CCO"))
    result, why = sp.template_with_hydrogens(residue_mol(real_template), other)
    assert result is None and "does not match" in why


def test_a_different_protonation_state_is_rejected(real_template):
    result, why = sp.template_with_hydrogens(residue_mol(real_template), ccd_for(real_template, extra_h=1))
    assert result is None and "hydrogens" in why


def test_entry_without_tables_is_rejected(real_template):
    empty = gemmi.cif.read_string("data_X\n_chem_comp.id X\n")[0]
    result, why = sp.template_with_hydrogens(residue_mol(real_template), empty)
    assert result is None and "no atom or bond table" in why


class FakeResponse:
    def __init__(self, text):
        self.text = text

    def read(self):
        return self.text.encode()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_dictionary_lookup_uses_the_expected_address(monkeypatch):
    seen = []

    def fake(url, timeout):
        seen.append((url, timeout))
        return FakeResponse("data_AQ4\n_chem_comp.id AQ4\n")

    monkeypatch.setattr(sp.urllib.request, "urlopen", fake)
    assert sp.fetch_ccd("aq4") is not None
    assert seen[0][0] == "https://files.rcsb.org/ligands/download/AQ4.cif" and seen[0][1] == 15


def test_dictionary_lookup_failures_mean_not_available(monkeypatch):
    def offline(url, timeout):
        raise OSError("network unreachable")

    monkeypatch.setattr(sp.urllib.request, "urlopen", offline)
    assert sp.fetch_ccd("AQ4") is None
    monkeypatch.setattr(sp.urllib.request, "urlopen", lambda url, timeout: FakeResponse("not a cif <<<"))
    assert sp.fetch_ccd("AQ4") is None


@pytest.mark.parametrize("name", ["../etc/passwd", "AQ4/../x", "", "A B", "TOOLONGNAME", "AQ4?x=1"])
def test_untrusted_residue_names_never_reach_the_network(monkeypatch, name):
    monkeypatch.setattr(sp.urllib.request, "urlopen", lambda *a, **k: pytest.fail("no request may be made"))
    assert sp.fetch_ccd(name) is None


# --------------------------------------------------------------------------
# Combined result files
# --------------------------------------------------------------------------


def docked_pose(energy, rank=1, shift=0.0, clusters=20):
    mol = Chem.AddHs(Chem.SDMolSupplier(LIGAND)[0], addCoords=True)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    mol.SetProp("meeko", json.dumps({"is_sidechain": [False], "free_energy": energy,
                                     "intermolecular_energy": energy - 1.0, "internal_energy": -1.0,
                                     "cluster_size": clusters, "cluster_id": rank, "rank_in_cluster": 1}))
    return mol


@pytest.fixture
def run_folder(tmp_path):
    """A finished run: one clean ligand with two cluster leads, one flagged ligand, one failure."""
    out = tmp_path / "out"
    for index, poses in ((0, [-9.7, -8.1]), (1, [-7.2])):
        folder = out / str(index)
        folder.mkdir(parents=True)
        with Chem.SDWriter(str(folder / "ligand_docked.sdf")) as w:
            for rank, e in enumerate(poses, 1):
                w.write(docked_pose(e, rank))
    (out / "2").mkdir()
    ad.write_summary(str(out), {
        0: {"index": 0, "name": "good ligand", "status": "docked", "anchor_max_dev": 0.8, "core_rmsd": 0.9, "note": ""},
        1: {"index": 1, "name": "odd/ligand", "status": "flagged", "anchor_max_dev": 4.1, "core_rmsd": 3.2,
            "note": "an anchor atom is 4.1 A from its target"},
        2: {"index": 2, "name": "broken", "status": "failed", "note": "mk_prepare_ligand.py failed (exit code 1)"}})
    return out


def test_scores_are_read_from_the_docked_pose(run_folder):
    rows = {r["index"]: r for r in results.collect(str(run_folder))}
    assert rows["0"]["free_energy"] == -9.7 and rows["0"]["n_poses"] == 2 and rows["0"]["cluster_size"] == 20
    assert rows["0"]["intermolecular_energy"] == pytest.approx(-10.7)
    assert rows["2"]["status"] == "failed" and rows["2"]["n_poses"] == 0 and rows["2"]["free_energy"] == ""


def test_results_table_is_a_csv_with_one_row_per_ligand(run_folder):
    results.write_all(str(run_folder))
    with open(run_folder / "results.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert [r["name"] for r in rows] == ["good ligand", "odd/ligand", "broken"]
    assert rows[0]["free_energy"] == "-9.7" and rows[0]["anchor_max_dev"] == "0.8"
    assert rows[2]["free_energy"] == "" and "mk_prepare_ligand" in rows[2]["note"]
    assert list(rows[0]) == list(results.CSV_COLUMNS)


def test_top_poses_hold_one_complete_molecule_per_docked_ligand(run_folder):
    made = results.write_all(str(run_folder))
    top = [m for m in Chem.SDMolSupplier(str(run_folder / "top_poses.sdf"), removeHs=False) if m]
    assert made["docked"] == 2 and len(top) == 2  # the failed ligand has no pose
    first = top[0]
    assert first.GetProp("_Name") == "good ligand" and first.GetProp("free_energy") == "-9.7"
    assert first.GetProp("pose_rank") == "1" and first.GetProp("status") == "docked"
    assert first.GetProp("anchor_max_dev") == "0.8"
    assert sum(a.GetAtomicNum() == 1 for a in first.GetAtoms()) == 23  # every hydrogen is there
    assert top[1].GetProp("status") == "flagged"


def test_all_poses_hold_every_cluster_lead(run_folder):
    made = results.write_all(str(run_folder))
    everything = [m for m in Chem.SDMolSupplier(str(run_folder / "all_poses.sdf"), removeHs=False) if m]
    assert made["poses"] == 3 and len(everything) == 3
    assert [m.GetProp("pose_rank") for m in everything] == ["1", "2", "1"]
    assert everything[1].GetProp("free_energy") == "-8.1" and not everything[1].HasProp("status")


def test_view_script_loads_everything_and_marks_flagged_poses(run_folder):
    results.write_all(str(run_folder))
    text = (run_folder / "view_all.pml").read_text()
    assert "load receptor_with_hydrogens.pdb, receptor" in text and "load template_with_hydrogens.sdf, template" in text
    assert "load 0/ligand_docked.sdf, L0_good_ligand" in text
    assert "load 1/ligand_docked.sdf, L1_odd_ligand" in text          # names are made safe for PyMOL
    assert "broken" not in text                                        # nothing to show for a failure
    assert "group ligands, L0_good_ligand L1_odd_ligand" in text
    assert "color orange, (L1_odd_ligand) and elem C" in text
    assert "save docking_all.pse" in text
    assert not any(l.rstrip().endswith(tuple("abcdefghijklmnopqrstuvwxyz")) and "#" in l and not l.startswith("#")
                   for l in text.splitlines())                        # no trailing comments after a command


def test_view_script_without_results_still_opens_the_template(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    ad.write_summary(str(out), {0: {"index": 0, "name": "x", "status": "failed", "note": "boom"}})
    results.write_all(str(out))
    text = (out / "view_all.pml").read_text()
    assert "zoom template" in text and "group ligands" not in text


def test_standalone_command(run_folder, capsys):
    assert results.main([str(run_folder)]) == 0
    assert "2 of 3 ligand(s) docked" in capsys.readouterr().out
    assert results.main([]) == 2 and results.main([str(run_folder / "nowhere")]) == 2


# --------------------------------------------------------------------------
# The docking script puts it together
# --------------------------------------------------------------------------


class Args:
    def __init__(self, tmp_path, receptor=RECEPTOR, template=TEMPLATE, view=None, add_hydrogens=False):
        self.docking_dir = str(tmp_path)
        self.receptor_file, self.template_file = receptor, template
        self.template_view_file, self.add_hydrogens = view, add_hydrogens


def test_docking_script_writes_the_display_files(tmp_path, capsys):
    ad.prepare_view_files(Args(tmp_path), "/adfr")
    assert results.count_hydrogens(str(tmp_path / results.VIEW_RECEPTOR)) == 2975
    template = Chem.SDMolSupplier(str(tmp_path / results.VIEW_TEMPLATE), removeHs=False)[0]
    assert sum(a.GetAtomicNum() == 1 for a in template.GetAtoms()) == 19
    assert "Receptor for viewing: 2975 hydrogens" in capsys.readouterr().out


def test_a_heavy_atom_template_from_a_structure_is_never_given_guessed_hydrogens(tmp_path):
    """Its bond orders are unknown, so RDKit would add the wrong hydrogens; the display copy is used instead."""
    from rdkit.Chem import rdmolops

    single = Chem.RWMol(Chem.SDMolSupplier(TEMPLATE)[0])
    for bond in single.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in single.GetAtoms():
        atom.SetIsAromatic(False)
    heavy = tmp_path / "heavy.sdf"
    heavy.write_text(Chem.MolToMolBlock(single.GetMol()) + "$$$$\n")
    guessed = tmp_path / "guessed.sdf"
    n_guess, _ = results.make_view_template(str(heavy), str(guessed))
    assert n_guess != 19  # guessing from single bonds is wrong, which is why the CLI always passes the display file
    view = tmp_path / "view.sdf"
    with Chem.SDWriter(str(view)) as w:
        w.write(Chem.AddHs(Chem.SDMolSupplier(TEMPLATE)[0], addCoords=True))
    n, _ = results.make_view_template(str(heavy), str(tmp_path / "t.sdf"), view_file=str(view))
    assert n == 19


def test_a_problem_with_the_display_files_does_not_stop_the_docking(tmp_path, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("ADFR crashed")

    monkeypatch.setattr(results, "make_view_receptor", boom)
    ad.prepare_view_files(Args(tmp_path), "/adfr")
    assert "Could not prepare the receptor with hydrogens" in capsys.readouterr().out
    assert (tmp_path / results.VIEW_RECEPTOR).is_file()  # the plain receptor is used instead


def test_per_ligand_script_uses_the_hydrogen_versions(tmp_path):
    ad.generate_pymol_session(str(tmp_path), TEMPLATE, "unused")
    text = (tmp_path / "view_results.pml").read_text()
    assert "load ../receptor_with_hydrogens.pdb, receptor" in text
    assert "load ../template_with_hydrogens.sdf, template" in text
    assert "load ligand_docked.sdf, docked_ligand" in text
