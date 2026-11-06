"""Tests for structure_prep.py and the wizard's 'start from a complex' path.

The complexes are built from the example Bace1 protein and template ligand, so no download is needed.
"""

import os
import subprocess
import sys

import gemmi
import numpy as np
import pytest
from rdkit import Chem, RDLogger

import anchored_docking as ad
import dock_cli as d
import structure_prep as sp

RDLogger.DisableLog("rdApp.*")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(REPO, "test", "Bace1_protein.pdb")
TEMPLATE = os.path.join(REPO, "test", "crystal_ligand.sdf")
LIGAND = os.path.join(REPO, "test", "CAT-13f.sdf")


@pytest.fixture(autouse=True)
def clear_warnings():
    d.WARNINGS.clear()
    yield
    d.WARNINGS.clear()


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    """The wizard looks the template up at rcsb.org; tests must never touch the network."""
    def refuse(name, timeout=15):
        return None

    monkeypatch.setattr(sp, "fetch_ccd", refuse)


def het_line(serial, name, resname, chain, resseq, xyz, element):
    x, y, z = xyz
    return (f"HETATM{serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00 20.00          {element:>2}\n")


def centre():
    return Chem.SDMolSupplier(TEMPLATE)[0].GetConformer().GetPositions().mean(axis=0)


def template_lines(chain="A", resseq=999, template_name="LIG"):
    mol = Chem.Mol(Chem.SDMolSupplier(TEMPLATE, removeHs=False)[0])
    lines = []
    for i, (atom, pos) in enumerate(zip(mol.GetAtoms(), mol.GetConformer().GetPositions())):
        label = f"{atom.GetSymbol()}{i + 1}"[:4]
        lines.append(het_line(8000 + i, f" {label:<3}", template_name, chain, resseq, pos, atom.GetSymbol()))
    return lines


def complex_lines(extra_zinc_far=True, second_chain=False, hydrogens=True, template_name="LIG"):
    c = centre()
    protein = [l for l in open(RECEPTOR) if l.startswith("ATOM")]
    if not hydrogens:
        protein = [l for l in protein if l[76:78].strip() != "H"]
    lines = list(protein) + ["TER\n"]
    if second_chain:  # a copy of the first chain far away, with its own copy of the ligand
        shifted = []
        for l in protein:
            x, y, z = float(l[30:38]) + 80, float(l[38:46]), float(l[46:54])
            shifted.append(f"{l[:21]}B{l[22:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
        lines += shifted + ["TER\n"]
    lines += template_lines(template_name=template_name)
    n = 9000
    rng = np.random.default_rng(0)
    for i in range(12):
        n += 1
        lines.append(het_line(n, " O  ", "HOH", "A", 2000 + i, c + rng.normal(scale=6, size=3), "O"))
    for k, (name, el, off) in enumerate([("S", "S", (0, 0, 0)), ("O1", "O", (1.4, 0, 0)), ("O2", "O", (-0.5, 1.3, 0)),
                                         ("O3", "O", (-0.5, -0.6, 1.2)), ("O4", "O", (-0.5, -0.6, -1.2))]):
        n += 1
        lines.append(het_line(n, f" {name:<3}", "SO4", "A", 3001, c + np.array([9, 0, 0]) + off, el))
    for k, (name, el, off) in enumerate([("C1", "C", (0, 0, 0)), ("O1", "O", (1.4, 0, 0)), ("C2", "C", (-0.7, 1.2, 0)),
                                         ("O2", "O", (-0.7, 2.5, 0.2)), ("C3", "C", (-1.4, -1, 0.5)),
                                         ("O3", "O", (-2.7, -1, 0.6))]):
        n += 1
        lines.append(het_line(n, f" {name:<3}", "GOL", "A", 3002, c + np.array([0, 10, 0]) + off, el))
    n += 1
    lines.append(het_line(n, "ZN  ", " ZN", "A", 3003, c + np.array([3.0, 3.0, 3.0]), "ZN"))  # near the pocket
    if extra_zinc_far:
        n += 1
        lines.append(het_line(n, "ZN  ", " ZN", "A", 3004, c + np.array([-40.0, 0, 0]), "ZN"))
    return lines + ["END\n"]


def write_complex(tmp_path, fmt="pdb", **kwargs):
    pdb = tmp_path / "complex.pdb"
    pdb.write_text("".join(complex_lines(**kwargs)))
    if fmt == "pdb":
        return str(pdb)
    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    cif = tmp_path / "complex.cif"
    st.make_mmcif_document().write_file(str(cif))
    return str(cif)


def load(path):
    st, n_models, altloc = sp.load_structure(path)
    return st, sp.inventory(st, n_models, altloc, path)


# --------------------------------------------------------------------------
# Reading and classifying
# --------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_inventory_classifies_everything_in_the_complex(tmp_path, fmt):
    st, inv = load(write_complex(tmp_path, fmt))
    assert list(inv.protein_chains) == ["A"] and inv.nucleic_chains == []
    assert inv.waters == 12
    kinds = {(h.name, h.kind) for h in inv.hets}
    assert {("LIG", "ligand"), ("SO4", "additive"), ("GOL", "additive"), ("ZN", "ion")} <= kinds
    assert inv.hydrogens > 0 and inv.protein_xyz.shape[1] == 3


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_the_bound_ligand_is_suggested_as_template(tmp_path, fmt):
    _, inv = load(write_complex(tmp_path, fmt))
    suggestion = sp.suggested_template(inv)
    assert suggestion.name == "LIG" and suggestion.heavy == 26
    candidates = sp.template_candidates(inv)
    assert candidates[0].name == "LIG"
    assert all(c.heavy >= sp.MIN_TEMPLATE_HEAVY for c in candidates)  # a lone zinc is never offered


def test_cofactors_are_not_suggested_over_a_real_ligand(tmp_path):
    lines = complex_lines()
    c = centre()
    for i in range(30):  # a large "NAD" touching the protein as much as anything (atom names must be unique)
        lines.insert(-1, het_line(7000 + i, f" C{i:<2}", "NAD", "A", 4000, c + np.array([2.0, 2.0, 2.0]) + 0.3 * i, "C"))
    path = tmp_path / "with_cofactor.pdb"
    path.write_text("".join(lines))
    _, inv = load(str(path))
    assert any(h.kind == "cofactor" for h in inv.hets)
    assert sp.suggested_template(inv).name == "LIG"


def test_unreadable_files_give_plain_errors(tmp_path):
    empty = tmp_path / "empty.pdb"
    empty.write_text("")
    with pytest.raises(sp.PrepError):
        sp.load_structure(str(empty))
    garbage = tmp_path / "garbage.cif"
    garbage.write_text("this is not a cif file\n")
    with pytest.raises(sp.PrepError):
        sp.load_structure(str(garbage))


def test_structure_without_protein_is_rejected(tmp_path):
    only_ligand = tmp_path / "ligand_only.pdb"
    only_ligand.write_text("".join(template_lines()) + "END\n")
    st, n_models, altloc = sp.load_structure(str(only_ligand))
    with pytest.raises(sp.PrepError, match="No protein"):
        sp.inventory(st, n_models, altloc)


def test_no_template_when_nothing_is_ligand_sized(tmp_path):
    protein = [l for l in open(RECEPTOR) if l.startswith("ATOM")]
    path = tmp_path / "apo.pdb"
    path.write_text("".join(protein) + het_line(9000, "ZN  ", " ZN", "A", 3000, (0, 0, 0), "ZN") + "END\n")
    _, inv = load(str(path))
    assert sp.template_candidates(inv) == [] and sp.suggested_template(inv) is None


# --------------------------------------------------------------------------
# Chains, molecules and the cleaned output
# --------------------------------------------------------------------------


def test_chains_near_the_template_are_the_default(tmp_path):
    st, inv = load(write_complex(tmp_path, second_chain=True))
    template = sp.suggested_template(inv)
    distances = sp.measure_distances(inv, template)
    assert distances["A"] < sp.CHAIN_NEAR < distances["B"]
    assert sp.default_chains(distances) == ["A"]
    assert sp.default_chains({"A": 30.0, "B": 25.0}) == ["B"]  # nothing near: the nearest chain


def test_plan_asks_only_about_relevant_molecules(tmp_path):
    st, inv = load(write_complex(tmp_path))
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    plan = sp.plan_molecules(inv, template, ["A"])
    assert [h.label for h in plan.ask] == ["ZN A3003"]  # the zinc by the pocket
    reasons = {h.label: why for h, why in plan.removed}
    assert reasons["SO4 A3001"] == "crystallisation additive" and reasons["GOL A3002"] == "crystallisation additive"
    assert "more than" in reasons["ZN A3004"]
    assert sp.default_keep(plan.ask[0]) is True  # within 5 A of the template


def test_other_copies_of_the_template_are_removed(tmp_path):
    st, inv = load(write_complex(tmp_path, second_chain=True))
    path = tmp_path / "two_ligands.pdb"
    lines = complex_lines(second_chain=True)
    shifted = [f"{l[:21]}B{l[22:30]}{float(l[30:38]) + 80:8.3f}{l[38:]}" for l in template_lines()]
    path.write_text("".join(lines[:-1] + shifted + ["END\n"]))
    st, inv = load(str(path))
    template = [h for h in inv.hets if h.name == "LIG" and h.chain == "A"][0]
    sp.measure_distances(inv, template)
    plan = sp.plan_molecules(inv, template, ["A", "B"])
    assert any(why == "another copy of the template ligand" for _, why in plan.removed)


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_cleaned_receptor_has_only_the_protein_and_what_was_kept(tmp_path, fmt):
    st, inv = load(write_complex(tmp_path, fmt))
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    plan = sp.plan_molecules(inv, template, ["A"])
    prepared = sp.prepare(st, inv, template, ["A"], [])
    het = {l[17:20].strip() for l in prepared.receptor_pdb.splitlines() if l.startswith("HETATM")}
    assert het == set()  # no template, waters, additives or ions
    assert sum(1 for l in prepared.receptor_pdb.splitlines() if l.startswith("ATOM")) == \
        sum(1 for l in open(RECEPTOR) if l.startswith("ATOM"))
    kept = sp.prepare(st, inv, template, ["A"], plan.ask)
    assert {l[17:20].strip() for l in kept.receptor_pdb.splitlines() if l.startswith("HETATM")} == {"ZN"}
    text = "\n".join(kept.summary)
    assert "Waters removed: 12" in text and "SO4" in text and "ZN A3004" in text and "LIG" in text


def test_second_chain_can_be_dropped_or_kept(tmp_path):
    st, inv = load(write_complex(tmp_path, second_chain=True))
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    one = sp.prepare(st, inv, template, ["A"], [])
    both = sp.prepare(st, inv, template, ["A", "B"], [])
    count = lambda p: sum(1 for l in p.receptor_pdb.splitlines() if l.startswith("ATOM"))
    assert count(both) == 2 * count(one)
    assert "Protein chains removed: B" in "\n".join(one.summary)


def test_alternate_locations_keep_one_conformation(tmp_path):
    lines = complex_lines()
    first = next(i for i, l in enumerate(lines) if l.startswith("ATOM"))
    original = lines[first]
    lines[first] = original[:16] + "A" + original[17:]
    lines.insert(first + 1, original[:16] + "B" + original[17:])
    path = tmp_path / "altloc.pdb"
    path.write_text("".join(lines))
    st, inv = load(str(path))
    assert inv.altloc_atoms == 2
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    prepared = sp.prepare(st, inv, template, ["A"], [])
    assert sum(1 for l in prepared.receptor_pdb.splitlines() if l.startswith("ATOM")) == \
        sum(1 for l in open(RECEPTOR) if l.startswith("ATOM"))
    assert "Alternate locations" in "\n".join(prepared.summary)


# --------------------------------------------------------------------------
# The extracted template
# --------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_template_keeps_the_coordinates_and_connectivity(tmp_path, fmt):
    st, inv = load(write_complex(tmp_path, fmt))
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    block = sp.prepare(st, inv, template, ["A"], []).template_molblock
    path = tmp_path / "template.sdf"
    path.write_text(block)
    mol = Chem.SDMolSupplier(str(path))[0]
    real = Chem.SDMolSupplier(TEMPLATE)[0]
    assert mol.GetNumHeavyAtoms() == real.GetNumHeavyAtoms() and mol.GetNumBonds() == real.GetNumBonds()
    assert mol.GetConformer().Is3D()
    a, b = mol.GetConformer().GetPositions(), real.GetConformer().GetPositions()
    assert np.abs(a - b).max() < 0.002  # only coordinate rounding (PDB keeps three decimals)


def test_extracted_template_gives_the_same_anchors_as_the_real_one(tmp_path):
    """Bond orders are lost on extraction; check that they do not matter for the shared core and anchors."""
    st, inv = load(write_complex(tmp_path))
    template = sp.suggested_template(inv)
    sp.measure_distances(inv, template)
    path = tmp_path / "template.sdf"
    path.write_text(sp.prepare(st, inv, template, ["A"], []).template_molblock)
    extracted = Chem.SDMolSupplier(str(path))[0]
    real = Chem.SDMolSupplier(TEMPLATE)[0]
    ligand = Chem.SDMolSupplier(LIGAND)[0]
    a = ad.find_anchor_points(real, ligand, 0, 100, None)
    b = ad.find_anchor_points(extracted, ligand, 0, 100, None)
    assert len(a.anchor_patterns) == len(b.anchor_patterns) == 17
    assert a.ptn2template == b.ptn2template
    assert np.abs(np.array(a.anchor_positions) - np.array(b.anchor_positions)).max() < 0.002


# --------------------------------------------------------------------------
# The wizard
# --------------------------------------------------------------------------


def save_choice():
    return "3" if subprocess.run(["which", "sbatch"], capture_output=True).returncode == 0 else "2"


def run_complex_wizard(monkeypatch, answers):
    import builtins

    monkeypatch.setattr(sys, "argv", ["dock_cli.py"])
    iterator = iter(answers)
    monkeypatch.setattr(builtins, "input", lambda prompt="": next(iterator))
    d.main()


def tail(output_dir):
    """Answers from 'how to run' to the end for a save-only run into `output_dir`."""
    return [save_choice(), str(output_dir), "", "", "", "", ""]


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_wizard_cleans_a_complex_and_saves_a_job(monkeypatch, tmp_path, capsys, fmt):
    structure = write_complex(tmp_path, fmt)
    out = tmp_path / "work" / "output"
    run_complex_wizard(monkeypatch, [
        "1", structure,          # start from a complex
        "",                      # template: the suggestion
        "yes",                   # keep the zinc by the pocket
        LIGAND, "", "",          # ligands, recommended settings, clean intermediates
        *tail(out)])
    text = capsys.readouterr().out
    assert "Waters removed: 12" in text and "crystallisation additive" in text and "Kept in the receptor: ZN" in text
    assert "Taken from the structure: LIG" in text
    job = tmp_path / "work" / "submit_dock_CAT-13f.sh"
    script = job.read_text()
    receptor = tmp_path / "work" / "dock_CAT-13f_receptor.pdb"
    template = tmp_path / "work" / "dock_CAT-13f_template.sdf"
    assert f"--receptor_file {receptor}" in script and f"--template_file {template}" in script
    assert receptor.is_file() and template.is_file()  # saved beside the results, not in a temporary folder
    assert {l[17:20].strip() for l in receptor.read_text().splitlines() if l.startswith("HETATM")} == {"ZN"}
    assert Chem.SDMolSupplier(str(template))[0].GetNumHeavyAtoms() == 26
    assert open(structure).read() == open(structure).read()  # the input is only read
    assert not out.exists()


def test_wizard_leaves_the_input_structure_unchanged(monkeypatch, tmp_path):
    structure = write_complex(tmp_path)
    before = open(structure).read()
    run_complex_wizard(monkeypatch, ["1", structure, "", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    assert open(structure).read() == before


def test_wizard_asks_which_chains_to_keep(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path, second_chain=True)
    run_complex_wizard(monkeypatch, [
        "1", structure, "",
        "A B",                   # keep both chains instead of the suggested one
        "no",                    # the zinc: remove
        LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "several protein chains" in text and "Nearest to template" in text
    receptor = (tmp_path / "w" / "dock_CAT-13f_receptor.pdb").read_text()
    chains = {l[21] for l in receptor.splitlines() if l.startswith("ATOM")}
    assert chains == {"A", "B"}


def test_wizard_rejects_an_unknown_chain_and_accepts_the_default(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path, second_chain=True)
    run_complex_wizard(monkeypatch, ["1", structure, "", "Z", "", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    assert "Choose from: A, B" in capsys.readouterr().out
    receptor = (tmp_path / "w" / "dock_CAT-13f_receptor.pdb").read_text()
    assert {l[21] for l in receptor.splitlines() if l.startswith("ATOM")} == {"A"}


def test_wizard_can_choose_another_file_when_no_template_exists(monkeypatch, tmp_path, capsys):
    apo = tmp_path / "apo.pdb"
    protein = [l for l in open(RECEPTOR) if l.startswith("ATOM")]
    apo.write_text("".join(protein) + "END\n")
    good = write_complex(tmp_path)
    run_complex_wizard(monkeypatch, ["1", str(apo), good, "", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    assert "No molecule in this file is big enough" in capsys.readouterr().out
    assert (tmp_path / "w" / "submit_dock_CAT-13f.sh").is_file()


def test_wizard_falls_back_to_adfr_hydrogens_when_pdb2pqr_is_declined(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path, hydrogens=False)
    run_complex_wizard(monkeypatch, [
        "1", structure, "", "no",
        "no",                    # protonate with PDB2PQR: no
        "",                      # add hydrogens with ADFR: yes (default)
        LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "no hydrogen atoms" in text and "quick rule-based guess" in text
    assert "--add_hydrogens" in (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()
    assert not (tmp_path / "w" / "dock_CAT-13f_protonation_report.txt").exists()


def test_wizard_can_decline_all_hydrogens_with_a_warning(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path, hydrogens=False)
    run_complex_wizard(monkeypatch, ["1", structure, "", "no", "no", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    assert "none will be added" in capsys.readouterr().out
    assert "--add_hydrogens" not in (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()


def test_hydrogens_are_not_offered_when_present(monkeypatch, tmp_path, capsys):
    run_complex_wizard(monkeypatch, ["1", write_complex(tmp_path), "", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "Protonate the protein" not in text and "Add hydrogens automatically" not in text
    assert "they are used as they are" in text
    assert "--add_hydrogens" not in (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()
    assert not (tmp_path / "w" / "dock_CAT-13f_protonation_report.txt").exists()


def test_wizard_protonates_with_pdb2pqr_by_default_at_ph_7_4(monkeypatch, tmp_path, capsys, fake_pdb2pqr):
    structure = write_complex(tmp_path, hydrogens=False)
    run_complex_wizard(monkeypatch, [
        "1", structure, "", "no",
        "",                      # protonate: yes (default)
        "",                      # pH: 7.4 (default)
        LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert fake_pdb2pqr.calls == [7.4]
    assert "Protonated at pH 7.4 with PDB2PQR and PROPKA" in text
    work = tmp_path / "w"
    receptor = (work / "dock_CAT-13f_receptor.pdb").read_text()
    assert " H1 " in receptor                                       # the protonated receptor is what is docked
    assert "--add_hydrogens" not in (work / "submit_dock_CAT-13f.sh").read_text()  # ADFR keeps these hydrogens
    report = (work / "dock_CAT-13f_protonation_report.txt").read_text()
    assert "pH: 7.4" in report and "Gaps in the chain" in report


def test_wizard_uses_the_ph_the_user_enters(monkeypatch, tmp_path, capsys, fake_pdb2pqr):
    structure = write_complex(tmp_path, hydrogens=False)
    run_complex_wizard(monkeypatch, ["1", structure, "", "no", "", "20", "6.5", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    assert fake_pdb2pqr.calls == [6.5]
    assert "pH: 6.5" in (tmp_path / "w" / "dock_CAT-13f_protonation_report.txt").read_text()


def test_wizard_falls_back_to_adfr_when_pdb2pqr_fails(monkeypatch, tmp_path, capsys, fake_pdb2pqr):
    fake_pdb2pqr.fail = True
    structure = write_complex(tmp_path, hydrogens=False)
    run_complex_wizard(monkeypatch, ["1", structure, "", "no", "", "", "", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "did not work" in text and "Falling back to ADFR" in text
    assert "--add_hydrogens" in (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()
    assert not (tmp_path / "w" / "dock_CAT-13f_protonation_report.txt").exists()


def test_wizard_shows_pocket_residues_that_deserve_a_look(monkeypatch, tmp_path, capsys, fake_pdb2pqr):
    structure = write_complex(tmp_path, hydrogens=False)
    # At pH 4 the catalytic aspartates (made-up pKa 3.8) are within the borderline window.
    run_complex_wizard(monkeypatch, ["1", structure, "", "no", "", "4", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "that deserve a look" in text and "ASP" in text and "3.8" in text
    assert "protonation_report.txt" in text


def test_separate_files_route_protonates_and_keeps_the_users_template(monkeypatch, tmp_path, capsys, fake_pdb2pqr):
    protein = tmp_path / "protein.pdb"
    protein.write_text("".join(l for l in open(RECEPTOR) if l.startswith("ATOM") and l[76:78].strip() != "H") + "END\n")
    run_complex_wizard(monkeypatch, ["2", str(protein), TEMPLATE, "", "", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    script = (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()
    assert f"--receptor_file {tmp_path / 'w' / 'dock_CAT-13f_receptor.pdb'}" in script
    assert f"--template_file {TEMPLATE}" in script
    assert "--template_view_file" not in script
    assert (tmp_path / "w" / "dock_CAT-13f_protonation_report.txt").is_file()
    assert " H1 " in (tmp_path / "w" / "dock_CAT-13f_receptor.pdb").read_text()
    assert protein.read_text().count(" H1 ") == 0                    # the user's own file is untouched


def test_gaps_next_to_the_template_are_warned_about(monkeypatch, tmp_path, capsys):
    """Cut residues out of the chain next to the pocket; the wizard must say so even with hydrogens present."""
    lines = complex_lines()
    atoms = [l for l in lines if l.startswith("ATOM")]
    c = centre()

    def dist(l):
        return np.linalg.norm(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) - c)
    near = min(atoms, key=dist)
    victim = int(near[22:26])
    kept = [l for l in lines if not (l.startswith("ATOM") and int(l[22:26]) in (victim, victim + 1))]
    path = tmp_path / "gapped.pdb"
    path.write_text("".join(kept))
    run_complex_wizard(monkeypatch, ["1", str(path), "", "no", LIGAND, "", "", *tail(tmp_path / "w" / "o")])
    text = capsys.readouterr().out
    assert "Chain gap between" in text and "from the template" in text
    assert any("Chain gap" in w for w in d.WARNINGS)


def test_output_folder_cannot_be_the_structure_folder(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path)
    run_complex_wizard(monkeypatch, [
        "1", structure, "", "no", LIGAND, "", "",
        save_choice(), str(tmp_path),           # contains the input structure: refused
        str(tmp_path / "w" / "o"), "", "", "", "", ""])
    assert "Choose another folder" in capsys.readouterr().out


def test_complex_flow_resumes_a_matching_earlier_run(monkeypatch, tmp_path, capsys):
    """The cleaned files are identical every time, so an earlier run with the same inputs can be resumed."""
    structure = write_complex(tmp_path)
    out = tmp_path / "w" / "output"
    answers = ["1", structure, "", "no", LIGAND, "", ""]
    run_complex_wizard(monkeypatch, answers + tail(out))
    saved_receptor = tmp_path / "w" / "dock_CAT-13f_receptor.pdb"
    saved_template = tmp_path / "w" / "dock_CAT-13f_template.sdf"
    out.mkdir()
    import json
    (out / "run_settings.json").write_text(json.dumps(
        ad.run_signature(str(saved_template), LIGAND, str(saved_receptor), d.max_anchor_count(), "spread", 20.0, 100,
                         None, None, None)))
    (out / "0").mkdir()
    (out / "0" / "ligand_docked.sdf").write_text("x")
    capsys.readouterr()
    run_complex_wizard(monkeypatch, answers + [save_choice(), str(out), "", "", "", "", ""])
    assert "same inputs and settings: 1 of 1" in capsys.readouterr().out
    assert "--resume" in (tmp_path / "w" / "submit_dock_CAT-13f.sh").read_text()


# --------------------------------------------------------------------------
# The docking script's hydrogen option
# --------------------------------------------------------------------------


def test_hydrogen_option_reaches_prepare_receptor(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(ad.subprocess, "run", lambda cmd, **kw: calls.append(cmd))
    ad.prepare_receptor(str(tmp_path), RECEPTOR, "/adfr")
    ad.prepare_receptor(str(tmp_path), RECEPTOR, "/adfr", add_hydrogens=True)
    assert "-A" not in calls[0]
    assert calls[1][-2:] == ["-A", "checkhydrogens"]


def test_hydrogen_option_is_part_of_the_resume_signature(tmp_path):
    base = ad.run_signature(TEMPLATE, LIGAND, RECEPTOR, 3, "random", 20.0, 100, None, None, None)
    with_h = ad.run_signature(TEMPLATE, LIGAND, RECEPTOR, 3, "random", 20.0, 100, None, None, None, True)
    assert ad.signature_differences(base, with_h) == ["add_hydrogens"]


def test_runs_saved_before_the_option_existed_still_resume(tmp_path):
    import json

    old = ad.run_signature(TEMPLATE, LIGAND, RECEPTOR, 3, "random", 20.0, 100, None, None, None)
    del old["add_hydrogens"]
    (tmp_path / "run_settings.json").write_text(json.dumps(old))
    now = ad.run_signature(TEMPLATE, LIGAND, RECEPTOR, 3, "random", 20.0, 100, None, None, None)
    assert ad.signature_differences(ad.read_run_settings(str(tmp_path)), now) == []
