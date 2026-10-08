"""Tests for the helper functions. Run them with:  pytest"""
import pytest

from brca_flow.tcga import (er_from_subtype, find_er_column, is_primary_tumor,
                            normalize_er_status, parse_os_status, patient_id,
                            sample_type_code)


def test_patient_id_from_long_barcode():
    assert patient_id("TCGA-A1-A0SP-01A-11R-A085-07") == "TCGA-A1-A0SP"


def test_patient_id_rejects_non_tcga():
    with pytest.raises(ValueError):
        patient_id("GTEX-1117F-0226")


def test_sample_type_codes():
    assert sample_type_code("TCGA-A1-A0SP-01A") == "01"
    assert sample_type_code("TCGA-A1-A0SP-11B") == "11"


def test_primary_tumor_detection():
    assert is_primary_tumor("TCGA-A1-A0SP-01A")
    assert not is_primary_tumor("TCGA-A1-A0SP-11A")


def test_sample_type_missing_field():
    with pytest.raises(ValueError):
        sample_type_code("TCGA-A1-A0SP")


@pytest.mark.parametrize("value, expected", [
    ("1:DECEASED", 1),
    ("0:LIVING", 0),
    ("  1:deceased ", 1),
    (None, None),
    ("[Not Available]", None),
    ("", None),
])
def test_parse_os_status(value, expected):
    assert parse_os_status(value) == expected


@pytest.mark.parametrize("value, expected", [
    ("Positive", "ER_positive"),
    ("negative", "ER_negative"),
    ("Indeterminate", None),
    ("[Not Evaluated]", None),
    (None, None),
])
def test_normalize_er_status(value, expected):
    assert normalize_er_status(value) == expected


def test_er_from_subtype():
    assert er_from_subtype("BRCA_LumA") == "ER_positive"
    assert er_from_subtype("BRCA_Basal") == "ER_negative"
    assert er_from_subtype("BRCA_Normal") is None


def test_find_er_column():
    assert find_er_column(["AGE", "ER_STATUS_BY_IHC"]) == "ER_STATUS_BY_IHC"
    assert find_er_column(["AGE", "SEX"]) is None
