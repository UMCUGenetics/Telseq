"""Data models for TelSeq."""

from dataclasses import dataclass, field
from typing import List

from telseq.parameters import GC_BIN_N, TEL_MOTIF_N


@dataclass
class ScanResults:
    """Results for a single read group."""

    sample: str = ""
    lib: str = ""
    bam: str = ""
    telcounts: List[int] = field(default_factory=lambda: [0] * TEL_MOTIF_N)
    gccounts: List[int] = field(default_factory=lambda: [0] * GC_BIN_N)
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
        self.telcounts = [0] * TEL_MOTIF_N
        self.gccounts = [0] * GC_BIN_N
        self.numTotal = 0
        self.numMapped = 0
        self.numDuplicates = 0
        self.telLenEstimate = 0.0
        self.n_exreadsExcluded = 0
        self.n_exreadsChrUnmatched = 0
        self.n_totalunfiltered = 0


@dataclass
class Range:
    """A genomic range."""

    start: int
    end: int


# Type aliases
StringVector = List[str]
RangeVector = List[Range]
