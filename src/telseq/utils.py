
"""Utility functions for TelSeq."""

import gzip
from typing import Dict, List, TextIO

from telseq.models import Range, RangeVector


def get_file_size(filename: str) -> int:
    """Return the size of a file in bytes."""
    try:
        with open(filename, "rb") as file:
            file.seek(0, 2)
            return file.tell()
    except OSError:
        try:
            with gzip.open(filename, "rb") as file:
                file.seek(0, 2)
                return file.tell()
        except OSError:
            return 0


def create_reader(filename: str) -> TextIO:
    """Create a file reader, handling gzip compression."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")

    return open(filename, "r")


def create_writer(filename: str) -> TextIO:
    """Create a file writer, handling gzip compression."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "wt")

    return open(filename, "w")


def read_bed_as_dict(filename: str) -> Dict[str, RangeVector]:
    """Read a BED file into a dictionary of chromosome intervals."""
    bed: Dict[str, RangeVector] = {}

    with create_reader(filename) as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) < 3:
                continue

            chromosome = parts[0]

            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            if start < 0 or end <= start:
                continue

            bed.setdefault(chromosome, []).append(
                Range(start=start, end=end)
            )

    return bed


def calc_gc(seq: str) -> float:
    """Calculate GC content of a sequence."""
    if not seq:
        return 0.0

    seq = seq.upper()
    num_gc = sum(1 for base in seq if base in "CG")

    return num_gc / len(seq)


def count_motif(
    seq: str,
    pattern: str,
    pattern_revcomp: str,
) -> int:
    """Count the maximum non-overlapping motif occurrences."""
    motif_count = count_occurrences(seq, pattern)
    reverse_count = count_occurrences(seq, pattern_revcomp)

    return max(motif_count, reverse_count)


def count_occurrences(seq: str, pattern: str) -> int:
    """Count non-overlapping occurrences of a pattern."""
    if not pattern:
        return 0

    count = 0
    position = 0
    pattern_length = len(pattern)

    while True:
        index = seq.find(pattern, position)

        if index == -1:
            break

        count += 1
        position = index + pattern_length

    return count


def ref_id_to_name(ref_id: int) -> str:
    """
    Convert a reference ID to a chromosome name.

    This assumes the reference order is chromosomes 1-22, X, Y.
    Prefer bam_file.get_reference_name() when scanning BAM files.
    """
    if 0 <= ref_id < 22:
        return str(ref_id + 1)

    if ref_id == 22:
        return "X"

    if ref_id == 23:
        return "Y"

    return "-1"
