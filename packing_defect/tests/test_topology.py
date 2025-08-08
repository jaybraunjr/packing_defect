from textwrap import dedent
from packing_defect.core.topology import TopologyReader
from packing_defect.core.classification import DefaultClassification
import json

def test_topology_reader_parses_resi_atom_lines(tmp_path):
    # minimal radii map and tiny topology block
    radii = {"CT2": 1.7, "OT": 1.5}
    top_contents = dedent("""
    * test
    RESI POPC 0.0
    ATOM C218 CT2  0.00
    ATOM O11  OT   0.00
    BOND C218 O11
    """).strip()
    topfile = tmp_path / "mini.rtf"
    topfile.write_text(top_contents)

    reader = TopologyReader(radii, DefaultClassification())
    out = reader.read("POPC", str(topfile))
    assert "C218" in out and "O11" in out
    r_c, code_c = out["C218"]
    r_o, code_o = out["O11"]
    assert r_c == 1.7 and r_o == 1.5
    # default classification: POPC C218 tail → 1, POPC O11 not tail → -1
    assert code_c == 1
    assert code_o == -1

def test_topology_reader_raises_on_missing_radius(tmp_path):
    radii = {"CT2": 1.7}
    topfile = tmp_path / "bad.rtf"
    topfile.write_text("RESI POPC\nATOM O11 OT 0.0\nBOND O11 O11\n")
    reader = TopologyReader(radii, DefaultClassification())
    try:
        reader.read("POPC", str(topfile))
        raised = False
    except KeyError:
        raised = True
    assert raised
