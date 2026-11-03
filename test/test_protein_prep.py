"""Tests for protein_prep.py: gap and missing-atom checks, log parsing, and protonation (fake and real PDB2PQR)."""

import os
import shutil

import numpy as np
import pytest
from rdkit import Chem

import protein_prep as pp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(REPO, "test", "Bace1_protein.pdb")
TEMPLATE = os.path.join(REPO, "test", "crystal_ligand.sdf")


def atom(serial, name, resname, chain, num, xyz, element=None):
    x, y, z = xyz
    element = element or name[0]
    return (f"ATOM  {serial:5d} {name:<4} {resname:>3} {chain}{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {element:>2}")


def tripeptide(chain="A", numbers=(1, 2, 3), gap_after_first=False, drop=()):
    """Three alanines in a line with 1.33 A peptide bonds (the middle one optionally far from the first)."""
    lines, serial = [], 1
    shift = 0.0
    for k, num in enumerate(numbers):
        if k == 1 and gap_after_first:
            shift = 5.0
        x = 3.8 * k + shift
        for name, dx, dy in (("N", 0.0, 0.0), ("CA", 1.2, 0.5), ("C", 2.4, 0.0), ("O", 2.6, -1.1), ("CB", 1.4, 1.9)):
            if (num, name) in drop:
                continue
            lines.append(atom(serial, name, "ALA", chain, num, (x + dx, dy, 0.0)))
            serial += 1
    return "\n".join(lines) + "\nEND\n"


# --------------------------------------------------------------------------
# Gaps and incomplete residues
# --------------------------------------------------------------------------


def test_a_continuous_chain_has_no_gaps():
    assert pp.find_gaps(pp.read_residues(tripeptide())) == []


def test_a_gap_is_found_and_measured():
    gaps = pp.find_gaps(pp.read_residues(tripeptide(gap_after_first=True)))
    assert len(gaps) == 1
    assert (gaps[0].before, gaps[0].after) == ("ALA A1", "ALA A2")
    assert gaps[0].bond_length > 5


def test_different_chains_are_not_a_gap():
    text = tripeptide(chain="A", numbers=(1, 2)) + tripeptide(chain="B", numbers=(1, 2))
    assert pp.find_gaps(pp.read_residues(text)) == []


def test_missing_side_chain_and_backbone_atoms_are_listed():
    found = pp.find_incomplete(pp.read_residues(tripeptide(drop={(2, "CB"), (3, "O")})))
    assert {(i.label, i.missing) for i in found} == {("ALA A2", ("CB",)), ("ALA A3", ("O",))}


def test_complete_residues_and_hetero_molecules_are_not_reported():
    text = tripeptide() + atom(99, "ZN", "ZN", "A", 100, (0, 0, 9), "ZN").replace("ATOM  ", "HETATM")
    assert pp.find_incomplete(pp.read_residues(text)) == []


def test_distance_to_the_template_is_reported_for_pocket_problems():
    template = np.array([[0.0, 0.0, 3.0]])
    found = pp.find_incomplete(pp.read_residues(tripeptide(drop={(1, "CB")})), template)
    assert found[0].distance == pytest.approx(3.0, abs=0.01)


def test_quality_lines_split_near_and_far_problems():
    residues = pp.read_residues(tripeptide(gap_after_first=True))
    near, far = np.array([[1.0, 0.0, 2.0]]), np.array([[100.0, 0.0, 0.0]])
    notes, urgent = pp.quality_lines(residues, near)
    assert urgent and not notes and "from the template" in urgent[0]
    notes, urgent = pp.quality_lines(residues, far)
    assert notes and not urgent and "should not matter" in notes[0]


def test_quality_lines_flag_incomplete_residues_next_to_the_template():
    residues = pp.read_residues(tripeptide(drop={(2, "CB")}))
    notes, urgent = pp.quality_lines(residues, np.array([[3.8 + 1.4, 1.9, 2.0]]))
    assert any("ALA A2" in u and "CB" in u for u in urgent)
    notes, urgent = pp.quality_lines(residues, np.array([[100.0, 0.0, 0.0]]))
    assert not urgent and "lack heavy atoms" in notes[0]


def test_clean_structures_give_no_lines():
    assert pp.quality_lines(pp.read_residues(tripeptide()), np.array([[0.0, 0.0, 2.0]])) == ([], [])


# --------------------------------------------------------------------------
# Log parsing (lines copied from a real PDB2PQR 3.8 run on 1M17)
# --------------------------------------------------------------------------

LOG = """\
INFO:Added atom OXT to residue PRO A 995
INFO:Added atom OXT to residue PRO A 995
WARNING:Unable to debump HIS A 811
WARNING:Missing atoms or failed protonation for LEU 977 A (BBN)
WARNING:Missing atoms or failed protonation for LEU 977 A (BBN)
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP 737 A     4.04       3.80
   ASP 872 A     7.59       3.80
   LYS 721 A    11.78      10.50
   N+  672 A     7.97       8.00
   C-  995 A     3.32       3.20
   CYS 773 B     9.90       9.00
--------------------------------------------------------------------------------------------------------
Free energy of   folding (kcal/mol) as a function of pH (using neutral reference)
   ASP 999 A     1.00       1.00
"""


def test_propka_values_are_read_from_the_summary_only():
    values = pp.parse_propka(LOG)
    assert values[("A", 872)] == ("ASP", 7.59)
    assert values[("A", 672)] == ("N+", 7.97) and values[("A", 995)] == ("C-", 3.32)
    assert values[("B", 773)] == ("CYS", 9.9)
    assert ("A", 999) not in values and len(values) == 6


def test_propka_with_no_summary_gives_nothing():
    assert pp.parse_propka("nothing useful\n") == {}


def test_repairs_and_problems_are_listed_once_each():
    repaired, problems = pp.parse_warnings(LOG)
    assert repaired == ["OXT added to PRO A995"]
    assert problems == ["HIS A811: hydrogens could not be placed without a clash", "LEU A977: protonation failed"]


# --------------------------------------------------------------------------
# protonate() with a fake PDB2PQR
# --------------------------------------------------------------------------


def with_hydrogen(text, rename=None, drop_atom=None, drop_residue=None):
    """What a runner returns: the input plus one hydrogen, optionally with a residue renamed or an atom lost."""
    lines = []
    for line in text.splitlines():
        if line.startswith("ATOM"):
            if drop_atom and (int(line[22:26]), line[12:16].strip()) == drop_atom:
                continue
            if drop_residue and int(line[22:26]) == drop_residue:
                continue
            if rename and int(line[22:26]) == rename[0]:
                line = line[:17] + rename[1] + line[20:]
        lines.append(line)
    first = next(l for l in lines if l.startswith("ATOM"))
    lines.insert(1, first[:12] + " H1 " + first[16:76] + " H" + first[78:])
    return "\n".join(lines) + "\n"


def runner(**kwargs):
    return lambda text, ph: (with_hydrogen(text, **kwargs), LOG)


def test_protonate_returns_the_new_structure_and_a_summary():
    prep = pp.protonate(tripeptide(), runner=runner(rename=(2, "HID")))
    assert prep.hydrogens == 1 and prep.ph == pp.DEFAULT_PH == 7.4
    assert "H1" in prep.pdb and prep.pdb.rstrip().endswith("END")
    assert "Protonated at pH 7.4" in prep.lines()[0]
    assert prep.repaired == ["OXT added to PRO A995"] and len(prep.problems) == 2


def test_residues_given_another_state_are_listed():
    asp = tripeptide().replace("ALA", "ASP")
    prep = pp.protonate(asp, runner=runner(rename=(2, "ASH")))
    assert prep.changed == [("ASP A2", "ASP", "ASH")]
    assert "ASH" in prep.report_text()


def test_heavy_atoms_that_moved_are_reported_as_flips():
    def mover(text, ph):
        out = []
        for line in with_hydrogen(text).splitlines():
            if line.startswith("ATOM") and int(line[22:26]) == 2 and line[12:16].strip() == "CB":
                line = line[:46] + f"{float(line[46:54]) + 2.0:8.3f}" + line[54:]
            out.append(line)
        return "\n".join(out) + "\n", LOG

    prep = pp.protonate(tripeptide(), runner=mover)
    assert prep.flipped == ["ALA A2"]
    assert "flipped" in " ".join(prep.lines()) and "ALA A2" in prep.report_text()
    assert pp.protonate(tripeptide(), runner=runner()).flipped == []     # small differences are not flips


def test_the_ph_is_passed_to_the_runner():
    seen = []

    def spy(text, ph):
        seen.append(ph)
        return with_hydrogen(text), LOG

    pp.protonate(tripeptide(), ph=6.0, runner=spy)
    assert seen == [6.0]


@pytest.mark.parametrize("ph", [-0.1, 14.5])
def test_ph_outside_the_scale_is_refused(ph):
    with pytest.raises(pp.ProteinPrepError, match="between 0 and 14"):
        pp.protonate(tripeptide(), ph=ph, runner=runner())


def test_a_file_without_protein_is_refused():
    with pytest.raises(pp.ProteinPrepError, match="no protein"):
        pp.protonate("HETATM    1 ZN    ZN A   1       0.000   0.000   0.000\nEND\n", runner=runner())


def test_output_that_lost_a_residue_is_rejected():
    with pytest.raises(pp.ProteinPrepError, match="number of residues"):
        pp.protonate(tripeptide(), runner=runner(drop_residue=3))


def test_output_that_lost_an_atom_is_rejected():
    with pytest.raises(pp.ProteinPrepError, match="dropped atoms.*CA"):
        pp.protonate(tripeptide(), runner=runner(drop_atom=(2, "CA")))


def test_output_without_hydrogens_is_rejected():
    with pytest.raises(pp.ProteinPrepError, match="no hydrogens"):
        pp.protonate(tripeptide(), runner=lambda text, ph: (text, LOG))


def test_header_records_are_dropped_before_running():
    seen = []

    def spy(text, ph):
        seen.append(text)
        return with_hydrogen(text), LOG

    pp.protonate("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n" + tripeptide(), runner=spy)
    assert not seen[0].startswith("CRYST1") and "CRYST1" not in seen[0]


def test_kept_molecules_without_hydrogens_are_reported():
    zinc = "HETATM  900 ZN    ZN A 100       0.000   0.000   9.000                       ZN\n"
    glycerol = ("HETATM  901  C1  GOL A 101       5.000   5.000   9.000\n"
                "HETATM  902  O1  GOL A 101       6.000   5.000   9.000\n")
    prep = pp.protonate(tripeptide().replace("END\n", "") + zinc + glycerol + "END\n", runner=runner())
    assert prep.kept_without_hydrogens == ["GOL"]       # a single-atom ion is not worth a warning
    assert "GOL" in prep.report_text()


def test_pocket_table_reports_titratable_residues_near_the_template():
    template = np.array([[0.0, 0.0, 2.0]])
    prep = pp.protonate(tripeptide(), template_xyz=template, runner=runner())
    assert prep.pocket == []                    # alanine does not titrate
    asp = tripeptide(numbers=(737, 738, 739)).replace("ALA", "ASP")
    prep = pp.protonate(asp, template_xyz=template, runner=runner())
    assert [p.label for p in prep.pocket] == ["ASP A737"]       # only residue 737 has a pKa in the log
    far = pp.protonate(asp, template_xyz=np.array([[100.0, 0.0, 0.0]]), runner=runner())
    assert far.pocket == []


def test_pocket_flags_residues_with_a_pka_near_the_ph():
    asp = tripeptide(numbers=(872, 873, 874)).replace("ALA", "ASP")
    prep = pp.protonate(asp, template_xyz=np.array([[0.0, 0.0, 2.0]]), runner=lambda t, p: (with_hydrogen(t), LOG))
    by_label = {p.label: p for p in prep.pocket}
    assert by_label["ASP A872"].pka == 7.59 and by_label["ASP A872"].notable        # 7.59 is close to pH 7.4
    assert "<-- look at this one" in prep.report_text()


def test_report_lists_gaps_and_incomplete_residues():
    prep = pp.protonate(tripeptide(gap_after_first=True, drop={(3, "CB")}), runner=runner())
    report = prep.report_text()
    assert "ALA A1 -> ALA A2" in report and "ALA A3: missing CB" in report


# --------------------------------------------------------------------------
# Real PDB2PQR (skipped where it is not installed)
# --------------------------------------------------------------------------

needs_pdb2pqr = pytest.mark.skipif(shutil.which("pdb2pqr") is None, reason="pdb2pqr is not installed")


@needs_pdb2pqr
def test_real_pdb2pqr_protonates_and_keeps_every_heavy_atom():
    text = "".join(l for l in open(RECEPTOR) if l.startswith("ATOM") and l[76:78].strip() != "H") + "END\n"
    template = Chem.SDMolSupplier(TEMPLATE)[0].GetConformer().GetPositions()
    prep = pp.protonate(text, template)
    before, after = pp.read_residues(text), pp.read_residues(prep.pdb)
    assert len(before) == len(after) and prep.hydrogens > 1000
    moved = set()
    for b, a in zip(before, after):
        heavy = [n for n in b.atoms if not pp.is_hydrogen(n)]
        assert set(heavy) <= set(a.atoms)
        if not np.allclose([b.atoms[n] for n in heavy], [a.atoms[n] for n in heavy], atol=0.01):
            moved.add(b.label)
    # Heavy atoms stay put except where PDB2PQR flipped an Asn, Gln or His, and it says which.
    assert moved <= set(prep.flipped) and all(label.split()[0] in ("ASN", "GLN", "HIS") for label in moved)
    assert set(prep.flipped) == moved
    assert prep.pocket and all(p.pka > 0 for p in prep.pocket)


@needs_pdb2pqr
def test_real_pdb2pqr_rebuilds_deleted_side_chain_atoms_and_reports_them():
    text = "".join(l for l in open(RECEPTOR) if l.startswith("ATOM") and l[76:78].strip() != "H") + "END\n"
    victim = next(l for l in text.splitlines() if l[12:16].strip() == "CD" and l[17:20] == "LYS")
    key = (int(victim[22:26]), victim[21])
    broken = "\n".join(l for l in text.splitlines()
                       if l.startswith("ATOM")
                       and not ((int(l[22:26]), l[21]) == key and l[12:16].strip() in ("CD", "CE", "NZ"))) + "\nEND\n"
    assert pp.find_incomplete(pp.read_residues(broken))
    prep = pp.protonate(broken)
    assert not [i for i in pp.find_incomplete(pp.read_residues(prep.pdb)) if i.label.endswith(str(key[0]))]
    assert any(f"{key[1]}{key[0]}" in r for r in prep.repaired)
    assert prep.incomplete                     # the input lacked atoms, and the report says so
