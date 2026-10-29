"""Utility functions for TelSeq."""

import gzip
from typing import Dict, List, Optional, TextIO, Tuple

from telseq.models import Range, RangeVector


def ref_id_to_name(ref_id: int) -> str:
    """Convert reference ID to chromosome name."""
    if ref_id < 22:
        return str(ref_id + 1)
    elif ref_id == 22:
        return "X"
    elif ref_id == 23:
        return "Y"
    else:
        return "-1"


def get_file_size(filename: str) -> int:
    """Return the size of a file in bytes."""
    try:
        with open(filename, "rb") as f:
            f.seek(0, 2)  # Seek to end
            return f.tell()
    except OSError:
        # Try gzip
        try:
            with gzip.open(filename, "rb") as f:
                f.seek(0, 2)
                return f.tell()
        except OSError:
            return 0


def create_reader(filename: str) -> TextIO:
    """Create a file reader, handling gzip compression."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")  # type: ignore[return-value]
    return open(filename, "r")


def create_writer(filename: str) -> TextIO:
    """Create a file writer, handling gzip compression."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "wt")  # type: ignore[return-value]
    return open(filename, "w")


def read_bed_as_dict(filename: str) -> Dict[str, RangeVector]:
    """Read a BED file and return a dictionary of chromosome to ranges."""
    bed: Dict[str, RangeVector] = {}
    with create_reader(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                chrm = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    if chrm not in bed:
                        bed[chrm] = []
                    bed[chrm].append(Range(start, end))
                except ValueError:
                    continue
    return bed


def calc_gc(seq: str) -> float:
    """Calculate GC content of a sequence."""
    if not seq:
        return 0.0
    num_gc = sum(1 for base in seq.upper() if base in "CG")
    return num_gc / len(seq)


def count_motif(seq: str, pattern: str, pattern_revcomp: str) -> int:
    """Count the maximum number of pattern or its reverse complement in a sequence."""
    motif_count1 = count_occurrences(seq, pattern)
    motif_count2 = count_occurrences(seq, pattern_revcomp)
    return max(motif_count1, motif_count2)


def count_occurrences(seq: str, pattern: str) -> int:
    """Count non-overlapping occurrences of pattern in seq."""
    count = 0
    pattern_len = len(pattern)
    if pattern_len == 0:
        return 0
    pos = 0
    while True:
        idx = seq.find(pattern, pos)
        if idx == -1:
            break
        count += 1
        pos = idx + pattern_len
    return count
