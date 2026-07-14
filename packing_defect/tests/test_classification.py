import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import json
from pathlib import Path

from packing_defect.core.classification import (
    DefaultClassification,
    UserDictClassification,
)


def test_default_cls():
    clf = DefaultClassification()

    # POPC tail atom → code 1; non-tail → -1
    assert clf.classify("POPC", "C22") == 1
    assert clf.classify("POPC", "P") == -1

    # TRIO glycerol atoms → 2; others → 3
    assert clf.classify("TRIO", "O11") == 2
    assert clf.classify("TRIO", "C99") == 3

    # Unknown residue → -1
    assert clf.classify("XXXX", "C22") == -1


def test_userdict_json(tmp_path: Path):
    data = {
        "POPC": {"C22": "tails", "P": "heads"},
        "TRIO": {"O11": "glycerol", "C31": "acyl"},
    }
    json_path = tmp_path / "rules.json"
    json_path.write_text(json.dumps(data))

    clf = UserDictClassification.from_json(str(json_path))

    # Consistent mapping for seen pairs
    c_tails = clf.classify("POPC", "C22")
    c_heads = clf.classify("POPC", "P")
    c_glyc = clf.classify("TRIO", "O11")
    c_acyl = clf.classify("TRIO", "C31")

    # All are positive nonzero codes
    assert all(x > 0 for x in (c_tails, c_heads, c_glyc, c_acyl))

    # Different labels map to different codes
    assert c_tails != c_heads
    assert c_glyc != c_acyl

    # Unlisted defaults to -1
    assert clf.classify("POPC", "C99") == -1
