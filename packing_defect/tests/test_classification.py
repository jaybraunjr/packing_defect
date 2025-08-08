from packing_defect.core.classification import DefaultClassification, UserDictClassification

def test_default_classification_popc_tail_and_other():
    dc = DefaultClassification()
    assert dc.classify("POPC", "C218") == 1
    assert dc.classify("POPC", "P") == -1  # headgroup

def test_default_classification_trio_glyc_vs_tail():
    dc = DefaultClassification()
    assert dc.classify("TRIO", "O11") == 2
    assert dc.classify("TRIO", "C216") in (2, 3)  # not glycerol set → 3

def test_userdictclassification_from_rules():
    rules = {("XYZ", "A1"): 42}
    uc = UserDictClassification(rules)
    assert uc.classify("XYZ", "A1") == 42
    assert uc.classify("XYZ", "A2") == -1
