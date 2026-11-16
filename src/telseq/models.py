
"""Data models for TelSeq."""

from dataclasses import dataclass, field
from typing import List

from telseq import parameters


@dataclass
class ScanResults:
    """Results for a single read group."""

    sample: str = ""
    lib: str = ""
    bam: str = ""

    telcounts: List[int] = field(
        default_factory=lambda: [0] * parameters.TEL_MOTIF_N
    )
    gccounts: List[int] = field(
        default_factory=lambda: [0] * parameters.GC_BIN_N
    )

    numTotal: int = 0
    numMapped: int = 0
    numDuplicates: int = 0
    telLenEstimate: float = 0.0

    # Logging data
    n_totalunfiltered: int = 0
    n_exreadsExcluded: int = 0
    n_exreadsChrUnmatched: int = 0

    def set_defaults(self) -> None:
        """Initialize counts to zeros."""
        self.telcounts = [0] * parameters.TEL_MOTIF_N
        self.gccounts = [0] * parameters.GC_BIN_N

        self.numTotal = 0
        self.numMapped = 0
        self.numDuplicates = 0
        self.telLenEstimate = 0.0

        self.n_totalunfiltered = 0
        self.n_exreadsExcluded = 0
        self.n_exreadsChrUnmatched = 0


@dataclass
class Range:
    """A genomic range using zero-based, half-open coordinates."""

    start: int
    end: int


# Type aliases
StringVector = List[str]
RangeVector = List[Range]
