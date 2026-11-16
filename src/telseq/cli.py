
"""Command-line interface for TelSeq."""

import argparse
import sys
from typing import Dict, List, Optional

from telseq import parameters
from telseq.models import Range, StringVector
from telseq.output import output_results, printlog
from telseq.scanner import Scanner
from telseq.utils import create_reader, get_file_size, read_bed_as_dict

PROGRAM_BIN = "telseq"
AUTHOR = "Zhihao Ding"

try:
    from telseq import __version__

    PROGRAM_VERSION = __version__
except ImportError:
    PROGRAM_VERSION = "0.0.1"

TELSEQ_VERSION_MESSAGE = (
    f"TelSeq Version {PROGRAM_VERSION}\n"
    f"Written by {AUTHOR}.\n"
    "Copyright 2013 Wellcome Trust Sanger Institute\n"
)

TELSEQ_USAGE_MESSAGE = (
    f"Program: TelSeq\n"
    f"Version: {PROGRAM_VERSION}\n"
    f"Contact: {AUTHOR} <zhihao.ding@gmail.com>\n\n"
    f"Usage: {PROGRAM_BIN} [OPTION] <in.1.bam> <in.2.bam> <...>\n"
    "Scan BAM and estimate telomere length.\n"
    "   <in.bam>                 one or more BAM files to be analysed. "
    "File names can also be passed from a pipe,\n"
    "                            with each row containing one BAM path.\n"
    "   -f, --bamlist=STR        a file that contains a list of file paths "
    "of BAMs. It should have only one column.\n"
    "                            -f has higher priority than <in.bam>. "
    "When specified, <in.bam> are ignored.\n"
    "   -o, --output=STR         output file for results. Ignored when "
    "input is from stdin, in which case output will be stdout.\n"
    "   -H                       remove header line, which is printed by default.\n"
    "   -h                       print the header line only.\n"
    "   -m                       merge read groups by taking a weighted "
    "average across read groups of a sample.\n"
    "   -u                       ignore read groups. Treat all reads in BAM "
    "as if they were from a same read group.\n"
    "   -k                       threshold of the amount of TTAGGG/CCCTAA "
    "repeats in read for a read to be considered telomeric. default = 7.\n"
    "   -r                       read length. default = 100\n"
    "   -z                       use user specified pattern for searching [ATGC]*.\n"
    "   -e, --exomebed=STR       specify exome regions in BED format. "
    "These regions will be excluded.\n"
    "   -w                       consider BAMs in the specified bamlist "
    "as one single BAM.\n"
)


class Headers:
    """Output headers."""

    def __init__(self) -> None:
        self.headers = self._generate_headers()

    @staticmethod
    def _generate_headers() -> List[str]:
        """Generate output headers using the current scan parameters."""
        headers = [
            parameters.LABEL_RG,
            parameters.LABEL_LB,
            parameters.LABEL_SAMPLE,
            parameters.LABEL_TOTAL,
            parameters.LABEL_MAPPED,
            parameters.LABEL_DUP,
            parameters.LABEL_LEN,
        ]

        for i in range(parameters.TEL_MOTIF_N):
            headers.append(f"{parameters.LABEL_TEL}{i}")

        for i in range(parameters.GC_BIN_N):
            headers.append(f"{parameters.LABEL_GC}{i}")

        return headers


def parse_options(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        prog=PROGRAM_BIN,
        description="Estimate telomere length from whole genome sequencing data (BAMs)",
        add_help=False,
    )

    parser.add_argument(
        "bam_files",
        nargs="*",
        help="One or more BAM files to be analysed",
    )

    parser.add_argument(
        "-f",
        "--bamlist",
        type=str,
        help="A file containing a list of BAM file paths",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Output file for results",
    )
    parser.add_argument(
        "-H",
        action="store_true",
        help="Remove header line",
    )
    parser.add_argument(
        "-h",
        action="store_true",
        help="Print header line only",
    )
    parser.add_argument(
        "-m",
        action="store_true",
        help="Merge read groups by weighted average",
    )
    parser.add_argument(
        "-u",
        action="store_true",
        help="Ignore read groups",
    )
    parser.add_argument(
        "-w",
        action="store_true",
        help="Consider BAMs in the specified bamlist as one single BAM",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=parameters.TEL_MOTIF_CUTOFF,
        help="Threshold for telomere motif repeats (default: 7)",
    )
    parser.add_argument(
        "-r",
        type=int,
        default=parameters.READ_LENGTH,
        help="Read length (default: 100)",
    )
    parser.add_argument(
        "-z",
        type=str,
        default="",
        help="User specified pattern for searching",
    )
    parser.add_argument(
        "-e",
        "--exomebed",
        type=str,
        default="",
        help="Exome regions in BED format to be excluded",
    )
    parser.add_argument(
        "--help",
        action="store_true",
        help="Display this help and exit",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Display version information and exit",
    )

    return parser.parse_args(args)


def process_bamlist_from_stdin() -> StringVector:
    """Read BAM paths from stdin."""
    bamlist: StringVector = []

    if not sys.stdin.isatty():
        for line in sys.stdin:
            line = line.strip()
            if line:
                bamlist.append(line)

    return bamlist


def process_bamlist_from_file(filename: str) -> StringVector:
    """Read BAM paths from a file."""
    filesize = get_file_size(filename)

    if filesize == 0:
        print(
            f"{PROGRAM_BIN}: BAMLIST file specified by -f is empty",
            file=sys.stderr,
        )
        sys.exit(1)

    bamlist: StringVector = []

    with create_reader(filename) as file_handle:
        for line in file_handle:
            line = line.strip()
            if line:
                bamlist.append(line)

    if not bamlist:
        print(
            f"{PROGRAM_BIN}: Could not find any sample in BAMLIST file specified.",
            file=sys.stderr,
        )
        sys.exit(1)

    return bamlist


def scan_bam() -> None:
    """Main scanning function."""
    options = parse_options()

    # Handle version and help before doing any scanning.
    if options.version:
        print(TELSEQ_VERSION_MESSAGE)
        return

    if options.help:
        print(TELSEQ_USAGE_MESSAGE)
        return

    # Update parameters before generating headers or creating the scanner.
    # This ensures that the motif bins and reverse complement are consistent.
    pattern = options.z if options.z else parameters.PATTERN

    try:
        parameters.update_parameters(
            read_length=options.r,
            pattern=pattern,
        )
    except ValueError as error:
        print(f"{PROGRAM_BIN}: {error}", file=sys.stderr)
        sys.exit(2)

    # Header-only mode must use the updated parameters.
    if options.h:
        headers = Headers()
        print(parameters.FIELD_SEP.join(headers.headers))
        return

    # Build BAM list. A BAM list file takes priority over positional inputs.
    if options.bamlist:
        bamlist = process_bamlist_from_file(options.bamlist)
    elif options.bam_files:
        bamlist = options.bam_files
    elif not sys.stdin.isatty():
        bamlist = process_bamlist_from_stdin()
    else:
        print(TELSEQ_USAGE_MESSAGE, file=sys.stderr)
        sys.exit(1)

    if not bamlist:
        print(TELSEQ_USAGE_MESSAGE, file=sys.stderr)
        sys.exit(1)

    # Load exome BED file if specified.
    exomebed: Dict[str, List[Range]] = {}

    if options.exomebed:
        exomebed = read_bed_as_dict(options.exomebed)
        print(
            f"loaded {len(exomebed)} exome regions\n",
            file=sys.stderr,
        )

    # Create scanner.
    scanner = Scanner(
        bam_files=bamlist,
        exome_bed=exomebed if exomebed else None,
        tel_k=options.k,
        read_length=parameters.READ_LENGTH,
        pattern=parameters.PATTERN,
        pattern_revcomp=parameters.PATTERN_REVCOMP,
        ignorerg=options.u,
        mergerg=options.m,
        onebam=options.w,
    )

    # Scan BAM files.
    resultlist = scanner.scan()

    # Write results.
    output_results(
        resultlist=resultlist,
        scanner=scanner,
        output_file=options.output,
        writerheader=not options.H,
        mergerg=options.m,
    )

    # Print exome filtering statistics when applicable.
    if exomebed:
        printlog(resultlist)

    print("Completed writing results", file=sys.stderr)


def main() -> None:
    """Main entry point."""
    scan_bam()


if __name__ == "__main__":
    main()
