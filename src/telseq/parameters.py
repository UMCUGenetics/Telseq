"""Scan parameters and constants for TelSeq."""

from typing import Final

# Field separator for output
FIELD_SEP: Final[str] = "\t"

# Genome length at telomeric GC (in base pairs)
GENOME_LENGTH_AT_TEL_GC: Final[int] = 332720800

# Output unit (1000 for kb, 1 for bp)
LENGTH_UNIT: Final[int] = 1000

# Default read length
READ_LENGTH: int = 100

# Telomere motif pattern and its reverse complement
PATTERN: str = "TTAGGG"
PATTERN_REVCOMP: str = "CCCTAA"

# Number of telomere ends (human: 46 chromosomes)
TELOMERE_ENDS: Final[int] = 46

# GC content bounds
GC_LOWERBOUND: Final[float] = 0.4
GC_UPPERBOUND: Final[float] = 0.6
GC_BINSIZE: Final[float] = 0.02

# Telomeric GC content bounds (for telomere length estimation)
GC_TELOMERIC_LOWERBOUND: Final[float] = 0.48
GC_TELOMERIC_UPPERBOUND: Final[float] = 0.52

# Maximum motif counts - calculated dynamically based on read length
# Add 1 to include 0 count
TEL_MOTIF_CUTOFF: Final[int] = 7

# Number of GC bins
GC_BIN_N: Final[int] = int((GC_UPPERBOUND - GC_LOWERBOUND) / GC_BINSIZE + 0.5)

# Maximum motif counts - will be calculated based on read length
TEL_MOTIF_N: int = READ_LENGTH // len(PATTERN) + 1

# Output labels
LABEL_RG: Final[str] = "ReadGroup"
LABEL_LB: Final[str] = "Library"
LABEL_SAMPLE: Final[str] = "Sample"
LABEL_TOTAL: Final[str] = "Total"
LABEL_MAPPED: Final[str] = "Mapped"
LABEL_DUP: Final[str] = "Duplicates"
LABEL_TEL: Final[str] = "TEL"
LABEL_GC: Final[str] = "GC"
LABEL_LEN: Final[str] = "LENGTH_ESTIMATE"


def update_parameters(read_length: int | None = None, pattern: str | None = None) -> None:
    """Update parameters based on read length and pattern."""
    global READ_LENGTH, PATTERN, PATTERN_REVCOMP, TEL_MOTIF_N

    if read_length is not None:
        READ_LENGTH = read_length
    if pattern is not None:
        PATTERN = pattern
        PATTERN_REVCOMP = reverse_complement(pattern)

    # Update total motif counts
    TEL_MOTIF_N = READ_LENGTH // len(PATTERN) + 1


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement_map = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
    }
    return "".join(complement_map.get(base, "N") for base in reversed(seq))
