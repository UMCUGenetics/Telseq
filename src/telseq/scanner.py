"""Core scanning logic for TelSeq."""

from typing import TYPE_CHECKING, Dict, List, Optional
import sys

if TYPE_CHECKING:
    import pysam
    HAS_PYSAM = True
else:
    pysam = None
    HAS_PYSAM = False

from telseq.models import Range, ScanResults, StringVector
from telseq.parameters import (
    GC_BINSIZE,
    GC_LOWERBOUND,
    GC_TELOMERIC_LOWERBOUND,
    GC_TELOMERIC_UPPERBOUND,
    GC_UPPERBOUND,
    GENOME_LENGTH_AT_TEL_GC,
    LENGTH_UNIT,
    PATTERN,
    PATTERN_REVCOMP,
    TEL_MOTIF_CUTOFF,
    TEL_MOTIF_N,
    TELOMERE_ENDS,
)
from telseq.utils import calc_gc, count_motif, ref_id_to_name


def _ensure_pysam() -> None:
    """Ensure pysam is imported and available."""
    global pysam, HAS_PYSAM
    if pysam is not None:
        return
    try:
        import pysam as pysam_module
        pysam = pysam_module
        HAS_PYSAM = True
    except ImportError:
        HAS_PYSAM = False
        raise ImportError(
            "pysam is required for BAM file reading. "
            "Install it with: uv pip install pysam"
        )


class Scanner:
    """Scanner for BAM files."""

    def __init__(
        self,
        bam_files: StringVector,
        exome_bed: Optional[Dict[str, List[Range]]] = None,
        tel_k: int = TEL_MOTIF_CUTOFF,
        read_length: int = 100,
        pattern: str = PATTERN,
        pattern_revcomp: str = PATTERN_REVCOMP,
        ignorerg: bool = False,
        writerheader: bool = True,
        mergerg: bool = False,
        onebam: bool = False,
    ) -> None:
        self.bam_files = bam_files
        self.exome_bed = exome_bed or {}
        self.tel_k = tel_k
        self.read_length = read_length
        self.pattern = pattern
        self.pattern_revcomp = pattern_revcomp
        self.ignorerg = ignorerg
        self.writerheader = writerheader
        self.mergerg = mergerg
        self.onebam = onebam

    def scan(self) -> List[Dict[str, ScanResults]]:
        """Scan BAM files and return results."""
        _ensure_pysam()

        resultlist: List[Dict[str, ScanResults]] = []
        is_exome = bool(self.exome_bed)

        print(f"{self.bam_files}\n{len(self.bam_files)} BAMs", file=sys.stderr)

        for bam_path in self.bam_files:
            resultmap: Dict[str, ScanResults] = {}
            lastfound: Dict[str, int] = {}

            print(f"Start analysing BAM {bam_path}", file=sys.stderr)

            # Open BAM file
            bam_file = pysam.AlignmentFile(bam_path, "rb")

            # Get header
            header = bam_file.header
            rggroups = "RG" in header

            if self.ignorerg:
                print("Treat all reads in BAM as if they were from a same sample", file=sys.stderr)
                results = ScanResults(sample="UNKNOWN")
                resultmap["UNKNOWN"] = results
            else:
                readgroups: Dict[str, str] = {}
                readlibs: Dict[str, str] = {}

                if rggroups:
                    for read_group in header["RG"]:
                        rg_id = read_group["ID"]
                        readgroups[rg_id] = read_group.get("SM", "UNKNOWN")
                        readlibs[rg_id] = read_group.get("LB", "UNKNOWN")

                    print(f"Specified BAM has {len(readgroups)} read groups", file=sys.stderr)

                    for rg_id, sample in readgroups.items():
                        results = ScanResults(sample=sample, lib=readlibs.get(rg_id, "UNKNOWN"))
                        resultmap[rg_id] = results
                else:
                    print("Warning: can't find RG tag in the BAM header", file=sys.stderr)
                    print("Warning: treat all reads in BAM as if they were from a same sample", file=sys.stderr)
                    results = ScanResults(sample="UNKNOWN", lib="UNKNOWN")
                    resultmap["UNKNOWN"] = results

            # Process reads
            nprocessed = 0  # Number of reads analyzed
            ntotal = 0  # Number of reads scanned

            for record in bam_file:
                ntotal += 1
                tag = "UNKNOWN"

                if rggroups:
                    # Skip reads without RG tag
                    if "RG" not in record.tags:
                        print(
                            f"can't find RG tag for read at position {record.reference_id}:{record.reference_start}",
                            file=sys.stderr,
                        )
                        print("skip this read", file=sys.stderr)
                        continue
                    # Get RG tag value
                    for rg_tag, rg_value in record.tags:
                        if rg_tag == "RG":
                            tag = rg_value
                            break

                # Skip reads with readgroup not defined in BAM header
                if tag not in resultmap:
                    print(
                        f"RG tag {tag} for read at position {record.reference_id}:{record.reference_start} "
                        "doesn't exist in BAM header.",
                        file=sys.stderr,
                    )
                    continue

                # For exome, exclude reads mapped to exome regions
                if is_exome and record.reference_id != -1:
                    chrm = ref_id_to_name(record.reference_id)
                    if chrm != "-1":
                        rg = Range(record.reference_start, record.reference_end)
                        if chrm in self.exome_bed:
                            ranges = self.exome_bed[chrm]
                            # Check for overlap
                            overlap_found = False
                            for exome_range in ranges:
                                if not (rg.end <= exome_range.start or rg.start >= exome_range.end):
                                    overlap_found = True
                                    break
                            if overlap_found:
                                resultmap[tag].n_exreadsExcluded += 1
                                continue
                        else:
                            resultmap[tag].n_exreadsChrUnmatched += 1

                resultmap[tag].numTotal += 1

                if record.is_mapped:
                    resultmap[tag].numMapped += 1

                if record.is_duplicate:
                    resultmap[tag].numDuplicates += 1

                # Get sequence
                seq = record.query_sequence
                if seq is None:
                    continue

                gc = calc_gc(seq)
                ptn_count = count_motif(seq, self.pattern, self.pattern_revcomp)

                # When the read length exceeds 100bp, number of patterns might exceed the boundary
                if ptn_count > TEL_MOTIF_N - 1:
                    continue

                resultmap[tag].telcounts[ptn_count] += 1

                if gc >= GC_LOWERBOUND and gc <= GC_UPPERBOUND:
                    # Get index for GC bin
                    idx = int((gc - GC_LOWERBOUND) / GC_BINSIZE)
                    if idx < 0 or idx >= len(resultmap[tag].gccounts):
                        print(
                            f"{nprocessed} GC:{gc} telcounts:{ptn_count} GC bin index out of bound: {idx}",
                            file=sys.stderr,
                        )
                        # Clamp to valid range
                        idx = max(0, min(idx, len(resultmap[tag].gccounts) - 1))
                    resultmap[tag].gccounts[idx] += 1

                nprocessed += 1

                if nprocessed % 10000000 == 0:
                    print(f"[scan] processed {nprocessed} reads\n", file=sys.stderr)

            bam_file.close()

            # Consider each BAM separately
            resultlist.append(resultmap)

            print(f"[scan] total reads in BAM scanned {ntotal}", file=sys.stderr)
            print("Completed scanning BAM", file=sys.stderr)

        if self.onebam:
            resultlist = self.merge_results_by_readgroup(resultlist)

        return resultlist

    def merge_results_by_readgroup(
        self, resultlist: List[Dict[str, ScanResults]]
    ) -> List[Dict[str, ScanResults]]:
        """Merge results in result list into one."""
        mergedresults: Dict[str, ScanResults] = {}

        for resultmap in resultlist:
            for rg, result in resultmap.items():
                if rg not in mergedresults:
                    mergedresults[rg] = result
                else:
                    self.add_results(mergedresults[rg], result)

        return [mergedresults]

    def add_results(self, x: ScanResults, y: ScanResults) -> None:
        """Combine counts in two ScanResults objects."""
        x.numTotal += y.numTotal
        x.numMapped += y.numMapped
        x.numDuplicates += y.numDuplicates
        x.n_exreadsExcluded += y.n_exreadsExcluded
        x.n_exreadsChrUnmatched += y.n_exreadsChrUnmatched
        x.n_totalunfiltered += y.n_totalunfiltered

        for j in range(min(len(x.telcounts), len(y.telcounts))):
            x.telcounts[j] += y.telcounts[j]
        for k in range(min(len(x.gccounts), len(y.gccounts))):
            x.gccounts[k] += y.gccounts[k]

    def calc_tel_length(self, results: ScanResults) -> float:
        """Calculate telomere length from scan results."""
        acc = 0
        for i in range(self.tel_k, len(results.telcounts)):
            acc += results.telcounts[i]

        gc_tel = 0
        for i in range(len(results.gccounts)):
            gc1 = GC_LOWERBOUND + GC_BINSIZE * i
            gc2 = GC_LOWERBOUND + GC_BINSIZE * (i + 1)
            if gc1 >= GC_TELOMERIC_LOWERBOUND and gc2 <= GC_TELOMERIC_UPPERBOUND:
                gc_tel += results.gccounts[i]

        if gc_tel == 0:
            return -1

        return (acc / gc_tel) * GENOME_LENGTH_AT_TEL_GC / LENGTH_UNIT / TELOMERE_ENDS
