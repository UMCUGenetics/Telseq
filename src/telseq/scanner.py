
"""Core scanning logic for TelSeq."""

import sys
from bisect import bisect_left
from typing import Any, Dict, List, Optional, Tuple

from telseq import parameters
from telseq.models import Range, ScanResults, StringVector
from telseq.utils import calc_gc, count_motif


PreparedExomeBed = Dict[str, Tuple[List[int], List[Range]]]


def _ensure_pysam() -> Any:
    """Import pysam only when BAM scanning is requested."""
    try:
        import pysam
    except ImportError as exc:
        raise ImportError(
            "pysam is required for BAM file reading. "
            "Install it with: uv pip install pysam"
        ) from exc

    return pysam


class Scanner:
    """Scan BAM files and calculate telomere-related read counts."""

    def __init__(
        self,
        bam_files: StringVector,
        exome_bed: Optional[Dict[str, List[Range]]] = None,
        tel_k: int = parameters.TEL_MOTIF_CUTOFF,
        read_length: Optional[int] = None,
        pattern: Optional[str] = None,
        pattern_revcomp: Optional[str] = None,
        ignorerg: bool = False,
        writerheader: bool = True,
        mergerg: bool = False,
        onebam: bool = False,
    ) -> None:
        self.bam_files = bam_files
        self.exome_bed = exome_bed or {}

        self.read_length = (
            read_length
            if read_length is not None
            else parameters.READ_LENGTH
        )
        self.pattern = (
            pattern
            if pattern is not None
            else parameters.PATTERN
        )
        self.pattern_revcomp = (
            pattern_revcomp
            if pattern_revcomp is not None
            else parameters.PATTERN_REVCOMP
        )

        if self.read_length <= 0:
            raise ValueError("Read length must be greater than zero.")

        if not self.pattern:
            raise ValueError("The motif pattern cannot be empty.")

        if len(self.pattern) != len(self.pattern_revcomp):
            raise ValueError(
                "The motif and reverse-complement pattern must have "
                "the same length."
            )

        self.pattern = self.pattern.upper()
        self.pattern_revcomp = self.pattern_revcomp.upper()

        self.tel_motif_n = (
            self.read_length // len(self.pattern) + 1
        )

        self.tel_k = tel_k

        if not 1 <= self.tel_k < self.tel_motif_n:
            raise ValueError(
                f"k must be between 1 and {self.tel_motif_n - 1} "
                f"for read length {self.read_length} and "
                f"motif length {len(self.pattern)}."
            )

        self.ignorerg = ignorerg
        self.onebam = onebam

        # Kept for compatibility with the CLI.
        # Output-specific behavior is handled in output.py.
        self.writerheader = writerheader
        self.mergerg = mergerg

        self._exome_bed = self._prepare_exome_bed(self.exome_bed)

    @staticmethod
    def _prepare_exome_bed(
        exome_bed: Dict[str, List[Range]],
    ) -> PreparedExomeBed:
        """Sort and merge exome intervals and cache their start positions."""
        prepared: PreparedExomeBed = {}

        for chromosome, ranges in exome_bed.items():
            sorted_ranges = sorted(
                ranges,
                key=lambda interval: (
                    interval.start,
                    interval.end,
                ),
            )

            merged: List[Range] = []

            for interval in sorted_ranges:
                if not merged or interval.start > merged[-1].end:
                    merged.append(
                        Range(
                            start=interval.start,
                            end=interval.end,
                        )
                    )
                else:
                    merged[-1].end = max(
                        merged[-1].end,
                        interval.end,
                    )

            starts = [interval.start for interval in merged]

            prepared[chromosome] = (starts, merged)

        return prepared

    @staticmethod
    def _new_result(
        sample: str,
        lib: str,
        bam: str,
        tel_motif_n: int,
    ) -> ScanResults:
        """Create a result with arrays sized for the active parameters."""
        result = ScanResults(
            sample=sample,
            lib=lib,
            bam=bam,
        )

        result.telcounts = [0] * tel_motif_n
        result.gccounts = [0] * parameters.GC_BIN_N

        return result

    def _read_overlaps_exome(
        self,
        chromosome: str,
        start: int,
        end: int,
    ) -> Optional[bool]:
        """
        Check whether a read overlaps an exome interval.

        Returns:
            True: read overlaps an exome interval.
            False: chromosome exists in BED, but no overlap was found.
            None: chromosome is absent from the BED.
        """
        chromosome_data = self._exome_bed.get(chromosome)

        if chromosome_data is None:
            return None

        starts, intervals = chromosome_data

        if not intervals:
            return False

        # Find the last interval starting before the read ends.
        index = bisect_left(starts, end) - 1

        if index < 0:
            return False

        interval = intervals[index]

        # BED and pysam use zero-based, half-open coordinates.
        return interval.end > start and interval.start < end

    def scan(self) -> List[Dict[str, ScanResults]]:
        """Scan all BAM files and return results per BAM/read group."""
        pysam = _ensure_pysam()

        resultlist: List[Dict[str, ScanResults]] = []
        is_exome = bool(self._exome_bed)
        nprocessed = 0

        print(self.bam_files)
        print(f"{len(self.bam_files)} BAMs")

        for bam_path in self.bam_files:
            resultmap: Dict[str, ScanResults] = {}
            ntotal = 0

            print(
                f"[scan] scanning BAM: {bam_path}",
                file=sys.stderr,
            )

            with pysam.AlignmentFile(bam_path, "rb") as bam_file:
                header = bam_file.header.to_dict()
                read_groups = header.get("RG", [])

                # Ignore RG information completely when ignorerg is enabled.
                rggroups = bool(read_groups) and not self.ignorerg

                if self.ignorerg:
                    resultmap["UNKNOWN"] = self._new_result(
                        sample="UNKNOWN",
                        lib="UNKNOWN",
                        bam=bam_path,
                        tel_motif_n=self.tel_motif_n,
                    )

                elif rggroups:
                    for read_group in read_groups:
                        rg_id = read_group.get("ID")

                        if not rg_id:
                            continue

                        sample = read_group.get("SM", "UNKNOWN")
                        library = read_group.get("LB", "UNKNOWN")

                        resultmap[rg_id] = self._new_result(
                            sample=sample,
                            lib=library,
                            bam=bam_path,
                            tel_motif_n=self.tel_motif_n,
                        )

                else:
                    resultmap["UNKNOWN"] = self._new_result(
                        sample="UNKNOWN",
                        lib="UNKNOWN",
                        bam=bam_path,
                        tel_motif_n=self.tel_motif_n,
                    )

                for record in bam_file:
                    ntotal += 1

                    # Determine the read group.
                    tag = "UNKNOWN"

                    if rggroups:
                        if not record.has_tag("RG"):
                            print(
                                f"[scan] read without RG tag in {bam_path}",
                                file=sys.stderr,
                            )
                            continue

                        tag = record.get_tag("RG")

                    if tag not in resultmap:
                        print(
                            f"[scan] read group {tag!r} not found "
                            f"in BAM header ({bam_path})",
                            file=sys.stderr,
                        )
                        continue

                    result = resultmap[tag]

                    # Exclude reads overlapping exome regions.
                    if is_exome and not record.is_unmapped:
                        chromosome = bam_file.get_reference_name(
                            record.reference_id
                        )
                        start = record.reference_start
                        end = record.reference_end

                        if (
                            chromosome is not None
                            and start >= 0
                            and end is not None
                        ):
                            overlap = self._read_overlaps_exome(
                                chromosome,
                                start,
                                end,
                            )

                            if overlap is True:
                                result.n_exreadsExcluded += 1
                                continue

                            if overlap is None:
                                result.n_exreadsChrUnmatched += 1

                    # Counters are incremented after exome filtering.
                    result.numTotal += 1

                    if not record.is_unmapped:
                        result.numMapped += 1

                    if record.is_duplicate:
                        result.numDuplicates += 1

                    sequence = record.query_sequence

                    if sequence is None:
                        continue

                    sequence = sequence.upper()

                    gc = calc_gc(sequence)

                    motif_count = count_motif(
                        sequence,
                        self.pattern,
                        self.pattern_revcomp,
                    )

                    # Ignore reads whose motif count exceeds the
                    # available motif-count bins.
                    if motif_count >= self.tel_motif_n:
                        continue

                    result.telcounts[motif_count] += 1

                    # Count GC content in the configured range.
                    if (
                        parameters.GC_LOWERBOUND
                        <= gc
                        <= parameters.GC_UPPERBOUND
                    ):
                        index = int(
                            (gc - parameters.GC_LOWERBOUND)
                            / parameters.GC_BINSIZE
                        )

                        # The upper GC boundary belongs to the final bin.
                        index = min(
                            index,
                            len(result.gccounts) - 1,
                        )

                        if index < 0:
                            raise ValueError(
                                f"Invalid GC bin index: {index}"
                            )

                        result.gccounts[index] += 1

                    nprocessed += 1

                    if nprocessed % 10_000_000 == 0:
                        print(
                            f"[scan] processed {nprocessed} reads",
                            file=sys.stderr,
                        )

            resultlist.append(resultmap)

            print(
                f"[scan] total reads in BAM scanned: {ntotal}",
                file=sys.stderr,
            )
            print(
                "[scan] completed scanning BAM",
                file=sys.stderr,
            )

        if self.onebam:
            resultlist = self.merge_results_by_readgroup(resultlist)

        return resultlist

    def merge_results_by_readgroup(
        self,
        resultlist: List[Dict[str, ScanResults]],
    ) -> List[Dict[str, ScanResults]]:
        """Merge results from multiple BAMs by read-group ID."""
        mergedresults: Dict[str, ScanResults] = {}

        for resultmap in resultlist:
            for rg, result in resultmap.items():
                if rg not in mergedresults:
                    mergedresults[rg] = result
                else:
                    self.add_results(
                        mergedresults[rg],
                        result,
                    )

        return [mergedresults]

    @staticmethod
    def add_results(
        x: ScanResults,
        y: ScanResults,
    ) -> None:
        """Add the counts from one ScanResults object to another."""
        x.numTotal += y.numTotal
        x.numMapped += y.numMapped
        x.numDuplicates += y.numDuplicates

        x.n_totalunfiltered += y.n_totalunfiltered
        x.n_exreadsExcluded += y.n_exreadsExcluded
        x.n_exreadsChrUnmatched += y.n_exreadsChrUnmatched

        if len(x.telcounts) != len(y.telcounts):
            raise ValueError(
                "Cannot merge results with different motif-bin sizes."
            )

        if len(x.gccounts) != len(y.gccounts):
            raise ValueError(
                "Cannot merge results with different GC-bin sizes."
            )

        for index in range(len(x.telcounts)):
            x.telcounts[index] += y.telcounts[index]

        for index in range(len(x.gccounts)):
            x.gccounts[index] += y.gccounts[index]

    def calc_tel_length(self, result: ScanResults) -> float:
        """Calculate the telomere length estimate for a read group."""
        telomeric_reads = sum(
            result.telcounts[self.tel_k :]
        )

        gc_telomeric_reads = 0

        for index, count in enumerate(result.gccounts):
            gc_start = (
                parameters.GC_LOWERBOUND
                + parameters.GC_BINSIZE * index
            )
            gc_end = (
                parameters.GC_LOWERBOUND
                + parameters.GC_BINSIZE * (index + 1)
            )

            if (
                gc_start >= parameters.GC_TELOMERIC_LOWERBOUND
                and gc_end <= parameters.GC_TELOMERIC_UPPERBOUND
            ):
                gc_telomeric_reads += count

        if gc_telomeric_reads == 0:
            return -1

        return (
            (telomeric_reads / gc_telomeric_reads)
            * parameters.GENOME_LENGTH_AT_TEL_GC
            / parameters.LENGTH_UNIT
            / parameters.TELOMERE_ENDS
        )
