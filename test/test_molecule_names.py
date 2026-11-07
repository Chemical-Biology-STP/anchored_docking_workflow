"""How the molecules of a structure are described to the user: names, formulas and SMILES from the PDB's chemical
component dictionary, with the formula of the atoms in the file as the fallback."""
import builtins

import gemmi
import pytest

import dock_cli as d
import structure_prep as sp
import test_structure_prep as tsp

AQ4 = """data_AQ4
_chem_comp.id AQ4
_chem_comp.name '[6,7-BIS(2-METHOXY-ETHOXY)QUINAZOLINE-4-YL]-(3-ETHYNYLPHENYL)AMINE'
_chem_comp.formula 'C22 H23 N3 O4'
_chem_comp.pdbx_synonyms ERLOTINIB
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
AQ4 SMILES OpenEye 1.0 'COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C'
AQ4 SMILES_CANONICAL CACTVS 3.341 'COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC'
"""
GOL = """data_GOL
_chem_comp.id GOL
_chem_comp.name GLYCEROL
_chem_comp.formula 'C3 H8 O3'
_chem_comp.pdbx_synonyms 'GLYCERIN; PROPANE-1,2,3-TRIOL'
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
GOL SMILES_CANONICAL CACTVS 3.341 'C(C(CO)O)O'
"""
ZN = """data_ZN
_chem_comp.id ZN
_chem_comp.name 'ZINC ION'
_chem_comp.formula Zn
_chem_comp.pdbx_synonyms ?
"""
LONG_STEREO = "[4R-(4ALPHA,5ALPHA,6BETA,7BETA)]-HEXAHYDRO-5,6-DIHYDROXY-1,3-BIS[2-NAPHTHYL-METHYL]-4,7-BIS(PHENYLMETHYL)-2H-1,3-DIAZEPIN-2-ONE"


def block(text):
    return gemmi.cif.read_string(text)[0]


@pytest.fixture(autouse=True)
def clear_warnings():
    d.WARNINGS.clear()
    yield
    d.WARNINGS.clear()


# --------------------------------------------------------------------------
# Choosing a name
# --------------------------------------------------------------------------


@pytest.mark.parametrize("systematic, synonyms, expected", [
    ("[6,7-BIS(2-METHOXY-ETHOXY)QUINAZOLINE-4-YL]-(3-ETHYNYLPHENYL)AMINE", "ERLOTINIB", "Erlotinib"),
    ("4-(4-METHYL-PIPERAZIN-1-YLMETHYL)-N-[4-METHYL-3-(4-PYRIDIN-3-YL-PYRIMIDIN-2-YLAMINO)-PHENYL]-BENZAMIDE",
     "STI-571; IMATINIB", "Imatinib"),                                    # a lab code with digits is skipped
    ("GLYCEROL", "GLYCERIN; PROPANE-1,2,3-TRIOL", "Glycerol"),            # a short systematic name is kept
    ("PROTOPORPHYRIN IX CONTAINING FE", "HEME", "Heme"),
    ("2-acetamido-2-deoxy-beta-D-glucopyranose", ";N-acetyl-beta-D-glucosamine; 2-acetamido-2-deoxy-D-glucose;",
     "N-acetyl-beta-D-glucosamine"),
    ("chloromethane", "Chloro Methyl group", "chloromethane"),            # mixed case is left alone
    ("ADENOSINE-5'-TRIPHOSPHATE", "", "Adenosine-5'-triphosphate"),
    ("SULFATE ION", "?", "Sulfate ion"),
    ("", "", ""),
])
def test_names_a_chemist_recognises(systematic, synonyms, expected):
    assert sp.friendly_name("XXX", systematic, synonyms) == expected


def test_a_long_systematic_name_without_a_common_one_loses_its_stereo_prefix_and_is_shortened():
    name = sp.friendly_name("XK2", LONG_STEREO, "")
    assert name.startswith("Hexahydro-5,6-dihydroxy") and name.endswith("...") and len(name) <= sp.NAME_MAX + 4
    assert "ALPHA" not in name.upper().replace("HEXAHYDRO", "")


# --------------------------------------------------------------------------
# Reading a dictionary entry
# --------------------------------------------------------------------------


def test_an_entry_gives_name_formula_and_canonical_smiles():
    info = sp.component_from_block("AQ4", block(AQ4))
    assert info.name == "Erlotinib" and info.formula == "C22H23N3O4"
    assert info.smiles == "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"          # the canonical one, not the other
    assert info.systematic.startswith("[6,7-BIS(2-METHOXY-ETHOXY)")


def test_an_entry_without_descriptors_or_synonyms_still_works():
    info = sp.component_from_block("GOL", block(GOL))
    assert (info.name, info.formula, info.smiles) == ("Glycerol", "C3H8O3", "C(C(CO)O)O")
    ion = sp.component_from_block("ZN", block(ZN))                       # no descriptor table and synonyms of "?"
    assert (ion.name, ion.formula, ion.smiles) == ("Zinc ion", "Zn", "")


def test_no_entry_gives_an_empty_description():
    assert sp.component_from_block("XYZ", None) == sp.Component("XYZ")


# --------------------------------------------------------------------------
# Looking molecules up
# --------------------------------------------------------------------------


def test_each_residue_name_is_fetched_once_whatever_the_case():
    asked = []
    components = sp.Components(lambda name, timeout=15: asked.append(name) or block(GOL))
    components.fetch("gol")
    components.fetch("GOL")
    components.prefetch(["Gol", "GOL"])
    assert asked == ["gol"]


def test_a_batch_is_fetched_with_a_short_timeout_and_counts_what_it_found():
    seen = []

    def fetch(name, timeout=15):
        seen.append((name, timeout))
        return block(AQ4) if name == "AQ4" else None

    components = sp.Components(fetch)
    assert components.prefetch(["AQ4", "GOL", "AQ4", "ZN"]) == 2          # both AQ4 entries count; the others were not found
    assert sorted(seen) == [("AQ4", sp.LOOKUP_SECONDS), ("GOL", sp.LOOKUP_SECONDS), ("ZN", sp.LOOKUP_SECONDS)]
    assert components.describe("AQ4").name == "Erlotinib" and components.describe("GOL") == sp.Component("GOL")


def test_nothing_found_is_not_an_error():
    components = sp.Components(lambda name, timeout=15: None)
    assert components.prefetch(["AQ4", "GOL"]) == 0
    assert components.describe("AQ4").name == ""


def test_the_download_function_is_looked_up_when_needed(monkeypatch):
    components = sp.Components()
    monkeypatch.setattr(sp, "fetch_ccd", lambda name, timeout=15: block(GOL))
    assert components.describe("GOL").name == "" and components.fetch("GOL") is not None
    assert components.describe("GOL").name == "Glycerol"


# --------------------------------------------------------------------------
# The formula of the atoms in the file
# --------------------------------------------------------------------------


def test_the_formula_of_the_atoms_in_the_file_is_known_without_any_lookup(tmp_path):
    st, inv = tsp.load(tsp.write_complex(tmp_path))
    by_name = {h.name: h.heavy_formula for h in inv.hets}
    assert by_name["SO4"] == "O4S" and by_name["GOL"] == "C3O3" and by_name["ZN"] == "Zn"
    mol = tsp.Chem.SDMolSupplier(tsp.TEMPLATE)[0]
    counts = {}
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    expected = "".join(f"{e}{n if n > 1 else ''}" for e, n in sorted(counts.items(), key=lambda kv: (kv[0] != "C", kv[0])))
    assert by_name["LIG"] == expected


def residue_with(*elements):
    res = gemmi.Residue()
    res.name = "UNL"
    for i, element in enumerate(elements):
        atom = gemmi.Atom()
        atom.name = f"X{i}"
        atom.element = gemmi.Element(element)
        res.add_atom(atom)
    return res


@pytest.mark.parametrize("elements, expected", [
    (["O", "C", "N", "C", "C"], "C3NO"),
    (["Br", "C", "C", "B", "Cl", "O"], "C2BBrClO"),          # carbon first, then alphabetical, although B and Br sort before C
    (["S", "O", "O", "O", "O"], "O4S"),
    (["H", "C", "O"], "CO"),                                 # hydrogens are not counted
])
def test_the_heavy_atom_formula_puts_carbon_first(elements, expected):
    assert sp._heavy_formula(residue_with(*elements)) == expected


# --------------------------------------------------------------------------
# The table the user sees
# --------------------------------------------------------------------------


def make_het(name, kind, formula="C2O", chain="A", num=1, icode=" "):
    return sp.Het(chain, name, num, icode, 5, kind, None, contacts=7, heavy_formula=formula)


def test_rows_show_the_name_the_complete_formula_and_the_position_in_separate_columns():
    components = sp.Components(lambda name, timeout=15: {"AQ4": block(AQ4), "GOL": block(GOL)}.get(name))
    shown = [make_het("AQ4", "ligand", "C22N3O4"), make_het("GOL", "additive", "C3O3")]
    components.prefetch(["AQ4", "GOL"])
    rows = d.molecule_rows(shown, components)
    assert rows[0] == (1, "Erlotinib", "C22H23N3O4", "AQ4", "A", "1", "ligand", 5, 7)    # code, chain and number apart
    assert rows[1][1:6] == ("Glycerol", "C3H8O3", "GOL", "A", "1")


def test_copies_of_the_same_molecule_are_told_apart_by_chain_and_number():
    components = sp.Components(lambda name, timeout=15: block(GOL))
    components.prefetch(["GOL"])
    shown = [make_het("GOL", "additive", chain="B", num=503), make_het("GOL", "additive", chain="A", num=503),
             make_het("GOL", "additive", chain="A", num=504)]
    rows = d.molecule_rows(shown, components)
    assert [(r[3], r[4], r[5]) for r in rows] == [("GOL", "B", "503"), ("GOL", "A", "503"), ("GOL", "A", "504")]
    assert len({(r[4], r[5]) for r in rows}) == 3


def test_an_insertion_code_stays_with_the_residue_number():
    assert make_het("JBG", "ligand", num=503, icode="A").number == "503A"
    assert make_het("JBG", "ligand", num=503, icode=" ").number == "503"
    assert make_het("JBG", "ligand", num=-2).number == "-2"                     # negative numbers occur in some files
    components = sp.Components(lambda name, timeout=15: None)
    assert d.molecule_rows([make_het("JBG", "ligand", num=503, icode="A")], components)[0][5] == "503A"


def test_a_file_without_chain_names_shows_a_dash_not_a_blank_cell():
    components = sp.Components(lambda name, timeout=15: None)
    assert d.molecule_rows([make_het("JBG", "ligand", chain=" ")], components)[0][4] == "-"
    assert d.molecule_rows([make_het("JBG", "ligand", chain="")], components)[0][4] == "-"


def test_without_a_lookup_the_rows_say_so_and_give_the_formula_of_the_atoms_present():
    components = sp.Components(lambda name, timeout=15: None)
    components.prefetch(["AQ4"])
    row = d.molecule_rows([make_het("AQ4", "ligand", "C22N3O4")], components)[0]
    assert row[1] == "(name unavailable)" and row[2] == "C22N3O4 (no H)"


def run_choose(monkeypatch, inv, fetch):
    monkeypatch.setattr(builtins, "input", lambda prompt="": "q")
    components = sp.Components(fetch)
    d.choose_template(inv, components)
    return components


def test_the_table_names_the_molecules_and_prints_smiles_for_ligand_like_ones_only(monkeypatch, tmp_path, capsys):
    st, inv = tsp.load(tsp.write_complex(tmp_path))
    entries = {"LIG": block(AQ4), "GOL": block(GOL), "SO4": block(GOL.replace("GOL", "SO4"))}
    run_choose(monkeypatch, inv, lambda name, timeout=15: entries.get(name))
    text = capsys.readouterr().out
    header = next(line for line in text.splitlines() if "Molecule" in line)
    assert header.split() == ["#", "Molecule", "Formula", "Code", "Chain", "Residue", "No.", "Kind", "Heavy", "atoms",
                              "Protein", "contacts"]
    assert "LIG A999" not in text                                                  # no longer one combined cell
    assert "Erlotinib" in text and "C22H23N3O4" in text
    import re
    assert re.search(r"Erlotinib\s+C22H23N3O4\s+LIG\s+A\s+999\s+ligand", text)    # name, formula, code, chain, number, kind
    assert "SMILES of the ligand-like molecules" in text
    assert "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC" in text                   # the ligand's SMILES
    smiles_lines = text.split("SMILES of the ligand-like molecules")[1]
    assert smiles_lines.count("\n    ") == 1                                      # glycerol and sulfate are not listed
    assert "Full chemical names:" in text and "QUINAZOLINE-4-YL" in text


def test_copies_of_a_molecule_share_one_line_under_all_their_row_numbers():
    components = sp.Components(lambda name, timeout=15: {"AQ4": block(AQ4), "GOL": block(GOL)}.get(name))
    components.prefetch(["AQ4", "GOL"])
    shown = [make_het("AQ4", "ligand", chain="B"), make_het("GOL", "additive"), make_het("AQ4", "ligand", chain="A"),
             make_het("AQ4", "ligand", chain="C")]
    smiles = d.per_molecule_notes(shown, components, lambda info, het: info.smiles if het.kind == "ligand" else "")
    assert smiles == [("1, 3, 4", "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC")]           # one line, not three
    names = d.per_molecule_notes(shown, components, lambda info, het: info.name)
    assert names == [("1, 3, 4", "Erlotinib"), ("2", "Glycerol")]                         # order follows the table


def test_different_molecules_with_the_same_text_are_not_merged():
    components = sp.Components(lambda name, timeout=15: block(GOL))
    components.prefetch(["GOL", "EDO"])
    shown = [make_het("GOL", "ligand"), make_het("EDO", "ligand")]
    both = d.per_molecule_notes(shown, components, lambda info, het: info.smiles)
    assert both == [("1", "C(C(CO)O)O"), ("2", "C(C(CO)O)O")]


def test_the_printed_table_lists_each_distinct_molecule_once(monkeypatch, tmp_path, capsys):
    two_copies = tmp_path / "copies"
    two_copies.mkdir()
    st, inv = tsp.load(tsp.write_complex(two_copies, second_chain=True))
    entries = {"LIG": block(AQ4), "GOL": block(GOL), "SO4": block(GOL.replace("GOL", "SO4"))}
    run_choose(monkeypatch, inv, lambda name, timeout=15: entries.get(name))
    text = capsys.readouterr().out
    smiles = text.split("SMILES of the ligand-like molecules")[1]
    assert smiles.count("COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC") == 1
    assert text.count("[6,7-BIS(2-METHOXY-ETHOXY)QUINAZOLINE-4-YL]") == 1


def test_the_full_name_is_not_repeated_when_the_column_already_shows_it(monkeypatch, tmp_path, capsys):
    st, inv = tsp.load(tsp.write_complex(tmp_path))
    entries = {"GOL": block(GOL), "SO4": block(GOL.replace("GOL", "SO4")), "LIG": block(GOL.replace("GOL", "LIG"))}
    run_choose(monkeypatch, inv, lambda name, timeout=15: entries.get(name))
    assert "Full chemical names:" not in capsys.readouterr().out                  # GLYCEROL is shown as Glycerol already


def test_offline_the_user_is_told_why_the_names_are_missing(monkeypatch, tmp_path, capsys):
    st, inv = tsp.load(tsp.write_complex(tmp_path))
    run_choose(monkeypatch, inv, lambda name, timeout=15: None)
    text = capsys.readouterr().out
    assert "(name unavailable)" in text and "(no H)" in text
    assert "could not be looked up at the PDB" in text
    assert "SMILES of the ligand-like molecules" not in text and "Full chemical names" not in text


def test_the_wizard_downloads_nothing_twice(monkeypatch, tmp_path):
    asked = []
    monkeypatch.setattr(sp, "fetch_ccd", lambda name, timeout=15: asked.append(name))
    structure = tsp.write_complex(tmp_path)
    tsp.run_complex_wizard(monkeypatch, ["1", structure, "", "no", tsp.LIGAND, "", "", *tsp.tail(tmp_path / "w" / "o")])
    assert sorted(asked) == ["GOL", "LIG", "SO4"] and "ZN" not in asked
