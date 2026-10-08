"""Helpers for TCGA identifiers and clinical fields.

These small functions are kept separate from the scripts so that they can be
tested on their own (see tests/test_tcga.py).
"""


def patient_id(barcode):
    """Return the patient part of a TCGA barcode.

    Example: 'TCGA-A1-A0SP-01A-11R-A085-07' -> 'TCGA-A1-A0SP'
    """
    parts = str(barcode).strip().split("-")
    if len(parts) < 3 or parts[0] != "TCGA":
        raise ValueError(f"Not a TCGA barcode: {barcode!r}")
    return "-".join(parts[:3])


def sample_type_code(barcode):
    """Return the two-digit sample type code of a TCGA barcode.

    '01' means primary tumor and '11' means solid tissue normal.
    Example: 'TCGA-A1-A0SP-01A' -> '01'
    """
    parts = str(barcode).strip().split("-")
    if len(parts) < 4 or len(parts[3]) < 2:
        raise ValueError(f"Barcode has no sample type field: {barcode!r}")
    return parts[3][:2]


def is_primary_tumor(barcode):
    """True if the barcode is a primary tumor sample (code 01)."""
    return sample_type_code(barcode) == "01"


def parse_os_status(value):
    """Convert a cBioPortal overall-survival status into 1, 0, or None.

    '1:DECEASED' -> 1 (the event happened), '0:LIVING' -> 0 (censored).
    Anything unrecognized, including missing values, becomes None.
    """
    if value is None:
        return None
    text = str(value).strip().upper()
    if text.startswith("1:") or text == "DECEASED":
        return 1
    if text.startswith("0:") or text == "LIVING":
        return 0
    return None


def normalize_er_status(value):
    """Turn a receptor-status word into 'ER_positive', 'ER_negative', or None."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if text == "positive":
        return "ER_positive"
    if text == "negative":
        return "ER_negative"
    return None


def er_from_subtype(subtype):
    """Approximate ER status from a PAM50 subtype (a fallback, not exact).

    Luminal A and B tumors are mostly ER positive; Basal and HER2-enriched
    tumors are mostly ER negative. Anything else returns None.
    """
    if subtype is None:
        return None
    text = str(subtype).upper()
    if "LUMA" in text or "LUMB" in text:
        return "ER_positive"
    if "BASAL" in text or "HER2" in text:
        return "ER_negative"
    return None


def find_er_column(columns):
    """Find the clinical column that holds ER status, or None if not found."""
    preferred = ["ER_STATUS_BY_IHC", "ER_STATUS", "ESTROGEN_RECEPTOR_STATUS"]
    for name in preferred:
        if name in columns:
            return name
    return None
