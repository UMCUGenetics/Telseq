"""Output formatting for TelSeq."""

import sys
from typing import TYPE_CHECKING, Dict, List, TextIO

from telseq.models import ScanResults
from telseq.parameters import (
    FIELD_SEP,
    GC_BIN_N,
    LABEL_DUP,
    LABEL_GC,
    LABEL_LEN,
    LABEL_LB,
    LABEL_MAPPED,
    LABEL_RG,
    LABEL_SAMPLE,
    LABEL_TEL,
    LABEL_TOTAL,
    TEL_MOTIF_N,
)

if TYPE_CHECKING:
    from telseq.scanner import Scanner


class Headers:
    """Output headers."""

    def __init__(self) -> None:
        self.headers = self._generate_headers()

    def _generate_headers(self) -> List[str]:
        """Generate output headers."""
        headers = [
            LABEL_RG,
            LABEL_LB,
            LABEL_SAMPLE,
            LABEL_TOTAL,
            LABEL_MAPPED,
            LABEL_DUP,
            LABEL_LEN,
        ]
        for i in range(TEL_MOTIF_N):
            headers.append(f"{LABEL_TEL}{i}")
        for i in range(GC_BIN_N):
            headers.append(f"{LABEL_GC}{i}")
        return headers


def _get_scanner() -> "Scanner":
    """Lazy import of Scanner."""
    from telseq.scanner import Scanner
    return Scanner([], tel_k=7)  # Default tel_k


def printout(rg: str, result: ScanResults, writer: TextIO) -> None:
    """Print a single result line."""
    writer.write(f"{rg}{FIELD_SEP}")
    writer.write(f"{result.lib}{FIELD_SEP}")
    writer.write(f"{result.sample}{FIELD_SEP}")
    writer.write(f"{result.numTotal}{FIELD_SEP}")
    writer.write(f"{result.numMapped}{FIELD_SEP}")
    writer.write(f"{result.numDuplicates}{FIELD_SEP}")

    # Calculate telomere length
    scanner = _get_scanner()
    result.telLenEstimate = scanner.calc_tel_length(result)

    if result.telLenEstimate == -1:
        print("Telomere length estimate unknown. No read was found with telomeric GC composition.", file=sys.stderr)
        writer.write(f"UNKNOWN{FIELD_SEP}")
    elif result.telLenEstimate > 1000000:
        print("Telomere length estimate unknown. Possibly due to not enough representation of genome.", file=sys.stderr)
        writer.write(f"UNKNOWN{FIELD_SEP}")
    elif result.telLenEstimate == 0:
        print(f"Telomere length estimate unknown. No read contains more than 7 telomere repeats.", file=sys.stderr)
        writer.write(f"UNKNOWN{FIELD_SEP}")
    else:
        writer.write(f"{result.telLenEstimate}{FIELD_SEP}")

    for count in result.telcounts:
        writer.write(f"{count}{FIELD_SEP}")
    for count in result.gccounts:
        writer.write(f"{count}{FIELD_SEP}")
    writer.write("\n")


def output_results(
    resultlist: List[Dict[str, ScanResults]],
    output_file: str = "",
    writerheader: bool = True,
    mergerg: bool = False,
) -> None:
    """Output results to file or stdout."""
    # Determine output stream
    if output_file:
        writer: TextIO = open(output_file, "w")
        tostdout = False
    else:
        writer = sys.stdout
        tostdout = True

    try:
        if writerheader:
            hd = Headers()
            writer.write(FIELD_SEP.join(hd.headers) + "\n")

        mergedrs: ScanResults
        grpnames: str = ""

        for i in range(len(resultlist)):
            resultmap = resultlist[i]

            # If merge read groups, take weighted average for all measures
            domg = mergerg and len(resultmap) > 1

            if domg:
                mergedrs = ScanResults()
                mergedrs.set_defaults()

            for rg, result in resultmap.items():
                if domg:
                    if len(grpnames) == 0:
                        grpnames = rg
                    else:
                        grpnames += "|" + rg

                    mergedrs.sample = result.sample
                    mergedrs.numTotal += result.numTotal
                    mergedrs.numMapped += result.numMapped * result.numTotal
                    mergedrs.numDuplicates += result.numDuplicates * result.numTotal

                    scanner = _get_scanner()
                    tel_len = scanner.calc_tel_length(result)
                    mergedrs.telLenEstimate += tel_len * result.numTotal

                    for j in range(min(len(mergedrs.telcounts), len(result.telcounts))):
                        mergedrs.telcounts[j] += result.telcounts[j] * result.numTotal
                    for k in range(min(len(mergedrs.gccounts), len(result.gccounts))):
                        mergedrs.gccounts[k] += result.gccounts[k] * result.numTotal
                    continue
                else:
                    printout(rg, result, writer)

            # In this case calculate weighted average
            if domg:
                if mergedrs.numTotal > 0:
                    mergedrs.numMapped /= mergedrs.numTotal
                    mergedrs.numDuplicates /= mergedrs.numTotal
                    mergedrs.telLenEstimate /= mergedrs.numTotal

                    for j in range(len(mergedrs.telcounts)):
                        mergedrs.telcounts[j] /= mergedrs.numTotal
                    for k in range(len(mergedrs.gccounts)):
                        mergedrs.gccounts[k] /= mergedrs.numTotal

                    mergedrs.numTotal /= len(resultmap)

                printout(grpnames, mergedrs, writer)

    finally:
        if not tostdout:
            writer.close()


def printlog(resultlist: List[Dict[str, ScanResults]]) -> None:
    """Print log information."""
    for i in range(len(resultlist)):
        rmap = resultlist[i]
        for rg, result in rmap.items():
            print(f"BAM:{rg}", file=sys.stderr)
            print(f"\tchr ID unmatched reads: {result.n_exreadsChrUnmatched}", file=sys.stderr)
            print(f"\texome reads excluded: {result.n_exreadsExcluded}", file=sys.stderr)
