
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

# Telomeric GC content bounds
GC_TELOMERIC_LOWERBOUND: Final[float] = 0.48
GC_TELOMERIC_UPPERBOUND: Final[float] = 0.52

# Minimum number of telomere repeats for a read to be considered telomeric
TEL_MOTIF_CUTOFF: Final[int] = 7

# Number of GC bins
GC_BIN_N: Final[int] = int(
    (GC_UPPERBOUND - GC_LOWERBOUND) / GC_BINSIZE + 0.5
)

# Number of motif-count bins, including bin zero
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


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement_map = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
    }

    seq = seq.upper()

    invalid_bases = set(seq) - set(complement_map)
    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases in motif: {sorted(invalid_bases)}"
        )

    return "".join(complement_map[base] for base in reversed(seq))


def update_parameters(
    read_length: int | None = None,
    pattern: str | None = None,
) -> None:
    """Update parameters based on read length and motif pattern."""
    global READ_LENGTH
    global PATTERN
    global PATTERN_REVCOMP
    global TEL_MOTIF_N

    if read_length is not None:
        if read_length <= 0:
            raise ValueError("Read length must be greater than zero.")
        READ_LENGTH = read_length

    if pattern is not None:
        pattern = pattern.upper()

        if not pattern:
            raise ValueError("The motif pattern cannot be empty.")

        PATTERN = pattern
        PATTERN_REVCOMP = reverse_complement(pattern)

    if not PATTERN:
        raise ValueError("The motif pattern cannot be empty.")

    TEL_MOTIF_N = READ_LENGTH // len(PATTERN) + 1
