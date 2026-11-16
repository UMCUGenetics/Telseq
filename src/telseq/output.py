
"""Output formatting for TelSeq."""

import sys
from typing import TYPE_CHECKING, Dict, List, TextIO

from telseq import parameters
from telseq.models import ScanResults
from telseq.utils import create_writer

if TYPE_CHECKING:
    from telseq.scanner import Scanner


class Headers:
    """Output headers."""

    def __init__(self) -> None:
        self.headers = self._generate_headers()

    @staticmethod
    def _generate_headers() -> List[str]:
        """Generate output headers using the active parameters."""
        headers = [
            parameters.LABEL_RG,
            parameters.LABEL_LB,
            parameters.LABEL_SAMPLE,
            parameters.LABEL_TOTAL,
            parameters.LABEL_MAPPED,
            parameters.LABEL_DUP,
            parameters.LABEL_LEN,
        ]

        for index in range(parameters.TEL_MOTIF_N):
            headers.append(f"{parameters.LABEL_TEL}{index}")

        for index in range(parameters.GC_BIN_N):
            headers.append(f"{parameters.LABEL_GC}{index}")

        return headers


def printout(
    rg: str,
    result: ScanResults,
    writer: TextIO,
    scanner: "Scanner",
) -> None:
    """Print a single result line."""
    writer.write(f"{rg}{parameters.FIELD_SEP}")
    writer.write(f"{result.lib}{parameters.FIELD_SEP}")
    writer.write(f"{result.sample}{parameters.FIELD_SEP}")
    writer.write(f"{result.numTotal}{parameters.FIELD_SEP}")
    writer.write(f"{result.numMapped}{parameters.FIELD_SEP}")
    writer.write(f"{result.numDuplicates}{parameters.FIELD_SEP}")

    result.telLenEstimate = scanner.calc_tel_length(result)

    if result.telLenEstimate == -1:
        print(
            "Telomere length estimate unknown. "
            "No read was found with telomeric GC composition.",
            file=sys.stderr,
        )
        writer.write(f"UNKNOWN{parameters.FIELD_SEP}")

    elif result.telLenEstimate > 1_000_000:
        print(
            "Telomere length estimate unknown. "
            "Possibly due to not enough representation of genome.",
            file=sys.stderr,
        )
        writer.write(f"UNKNOWN{parameters.FIELD_SEP}")

    elif result.telLenEstimate == 0:
        print(
            "Telomere length estimate unknown. "
            f"No read contains at least {scanner.tel_k} "
            "telomere repeats.",
            file=sys.stderr,
        )
        writer.write(f"UNKNOWN{parameters.FIELD_SEP}")

    else:
        writer.write(
            f"{result.telLenEstimate}{parameters.FIELD_SEP}"
        )

    for count in result.telcounts:
        writer.write(f"{count}{parameters.FIELD_SEP}")

    for count in result.gccounts:
        writer.write(f"{count}{parameters.FIELD_SEP}")

    writer.write("\n")


def _merge_read_groups(
    resultmap: Dict[str, ScanResults],
) -> tuple[str, ScanResults]:
    """
    Merge read groups using a weighted average.

    The weight for each read group is its numTotal value.
    This preserves the original C++ output behavior.
    """
    mergedrs = ScanResults()
    mergedrs.set_defaults()

    group_names = "|".join(resultmap)

    total_weight = sum(
        result.numTotal for result in resultmap.values()
    )

    if not resultmap:
        return group_names, mergedrs

    if total_weight == 0:
        mergedrs.sample = next(
            iter(resultmap.values())
        ).sample
        return group_names, mergedrs

    number_of_groups = len(resultmap)

    weighted_mapped = 0.0
    weighted_duplicates = 0.0

    weighted_telcounts = [
        0.0
    ] * len(mergedrs.telcounts)

    weighted_gccounts = [
        0.0
    ] * len(mergedrs.gccounts)

    for result in resultmap.values():
        weight = result.numTotal

        mergedrs.sample = result.sample

        weighted_mapped += result.numMapped * weight
        weighted_duplicates += result.numDuplicates * weight

        if len(result.telcounts) != len(weighted_telcounts):
            raise ValueError(
                "Cannot merge read groups with different motif-bin sizes."
            )

        if len(result.gccounts) != len(weighted_gccounts):
            raise ValueError(
                "Cannot merge read groups with different GC-bin sizes."
            )

        for index, count in enumerate(result.telcounts):
            weighted_telcounts[index] += count * weight

        for index, count in enumerate(result.gccounts):
            weighted_gccounts[index] += count * weight

    mergedrs.numTotal = total_weight / number_of_groups
    mergedrs.numMapped = weighted_mapped / total_weight
    mergedrs.numDuplicates = weighted_duplicates / total_weight

    mergedrs.telcounts = [
        count / total_weight
        for count in weighted_telcounts
    ]

    mergedrs.gccounts = [
        count / total_weight
        for count in weighted_gccounts
    ]

    return group_names, mergedrs


def output_results(
    resultlist: List[Dict[str, ScanResults]],
    scanner: "Scanner",
    output_file: str = "",
    writerheader: bool = True,
    mergerg: bool = False,
) -> None:
    """Output results to a file or stdout."""
    if output_file:
        writer = create_writer(output_file)
        close_writer = True
    else:
        writer = sys.stdout
        close_writer = False

    try:
        if writerheader:
            headers = Headers()
            writer.write(
                parameters.FIELD_SEP.join(headers.headers) + "\n"
            )

        for resultmap in resultlist:
            if not resultmap:
                continue

            # Merge read groups within this BAM only.
            if mergerg and len(resultmap) > 1:
                group_names, mergedrs = _merge_read_groups(resultmap)

                printout(
                    group_names,
                    mergedrs,
                    writer,
                    scanner,
                )

            else:
                for rg, result in resultmap.items():
                    printout(
                        rg,
                        result,
                        writer,
                        scanner,
                    )

    finally:
        if close_writer:
            writer.close()


def printlog(
    resultlist: List[Dict[str, ScanResults]],
) -> None:
    """Print exome filtering statistics."""
    for resultmap in resultlist:
        for rg, result in resultmap.items():
            print(f"BAM:{rg}", file=sys.stderr)
            print(
                f"\tchr ID unmatched reads: "
                f"{result.n_exreadsChrUnmatched}",
                file=sys.stderr,
            )
            print(
                f"\texome reads excluded: "
                f"{result.n_exreadsExcluded}",
                file=sys.stderr,
            )
