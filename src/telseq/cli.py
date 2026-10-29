"""Command-line interface for TelSeq."""

import argparse
import sys
from typing import Dict, List, Optional

from telseq.models import Range, ScanResults, StringVector
from telseq.output import output_results, printlog
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
    PATTERN,
    PATTERN_REVCOMP,
    READ_LENGTH,
    TEL_MOTIF_CUTOFF,
    TEL_MOTIF_N,
    update_parameters,
)
from telseq.scanner import Scanner
from telseq.utils import create_reader, get_file_size, read_bed_as_dict

PROGRAM_BIN = "telseq"
AUTHOR = "Zhihao Ding"

try:
    from telseq import __version__
    PROGRAM_VERSION = __version__
except ImportError:
    PROGRAM_VERSION = "0.0.1"

TELSEQ_VERSION_MESSAGE = f"TelSeq Version {PROGRAM_VERSION}\nWritten by {AUTHOR}.\nCopyright 2013 Wellcome Trust Sanger Institute\n"

TELSEQ_USAGE_MESSAGE = (
    f"Program: TelSeq\n"
    f"Version: {PROGRAM_VERSION}\n"
    f"Contact: {AUTHOR} <zhihao.ding@gmail.com>\n\n"
    f"Usage: {PROGRAM_BIN} [OPTION] <in.1.bam> <in.2.bam> <...> \n"
    f"Scan BAM and estimate telomere length. \n"
    f"   <in.bam>                 one or more BAM files to be analysed. File names can also be passed from a pipe,\n"
    f"                            with each row containing one BAM path.\n"
    f"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should has only one column,\n"
    f"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified,\n"
    f"                            <in.bam> are ignored.\n"
    f"   -o, --output=STR         output file for results. Ignored when input is from stdin, in which case output will be stdout.\n"
    f"   -H                       remove header line, which is printed by default.\n"
    f"   -h                       print the header line only.\n"
    f"   -m                       merge read groups by taking a weighted average across read groups of a sample.\n"
    f"   -u                       ignore read groups. Treat all reads in BAM as if they were from a same read group.\n"
    f"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
    f"   -r                       read length. default = 100\n"
    f"   -z                       use user specified pattern for searching [ATGC]*.\n"
    f"   -e, --exomebed=STR       specify exome regions in BED format. These regions will be excluded.\n"
    f"   -w                       consider BAMs in the specified bamlist as one single BAM.\n"
)


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


def parse_options(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command line options."""
    parser = argparse.ArgumentParser(
        prog=PROGRAM_BIN,
        description="Estimate telomere length from whole genome sequencing data (BAMs)",
        add_help=False,
    )

    # Positional arguments for BAM files
    parser.add_argument(
        "bam_files",
        nargs="*",
        help="One or more BAM files to be analysed",
    )

    # Optional arguments
    parser.add_argument(
        "-f", "--bamlist",
        type=str,
        help="A file containing a list of BAM file paths",
    )
    parser.add_argument(
        "-o", "--output",
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
        default=TEL_MOTIF_CUTOFF,
        help="Threshold for telomere motif repeats (default: 7)",
    )
    parser.add_argument(
        "-r",
        type=int,
        default=READ_LENGTH,
        help="Read length (default: 100)",
    )
    parser.add_argument(
        "-z",
        type=str,
        default="",
        help="User specified pattern for searching",
    )
    parser.add_argument(
        "-e", "--exomebed",
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

    # Parse arguments
    options = parser.parse_args(args)

    return options


def process_bamlist_from_stdin() -> StringVector:
    """Read BAM list from stdin."""
    bamlist: StringVector = []
    if not sys.stdin.isatty():
        for line in sys.stdin:
            line = line.strip()
            if line:
                bamlist.append(line)
    return bamlist


def process_bamlist_from_file(filename: str) -> StringVector:
    """Read BAM list from file."""
    bamlist: StringVector = []
    filesize = get_file_size(filename)
    if filesize == 0:
        print(f"{PROGRAM_BIN}: BAMLIST file specified by -f is empty", file=sys.stderr)
        sys.exit(1)

    with create_reader(filename) as f:
        for line in f:
            line = line.strip()
            if line:
                bamlist.append(line)

    if not bamlist:
        print(f"{PROGRAM_BIN}: Could not find any sample in BAMLIST file specified.", file=sys.stderr)
        sys.exit(1)

    return bamlist


def scan_bam() -> None:
    """Main scanning function."""
    options = parse_options()

    # Handle version and help
    if options.version:
        print(TELSEQ_VERSION_MESSAGE)
        sys.exit(0)
    if options.help:
        print(TELSEQ_USAGE_MESSAGE)
        sys.exit(0)

    # Handle header only
    if options.h:
        hd = Headers()
        print(FIELD_SEP.join(hd.headers))
        sys.exit(0)

    # Update parameters based on options
    if options.r != READ_LENGTH or options.z:
        pattern = options.z if options.z else PATTERN
        update_parameters(read_length=options.r, pattern=pattern)

    # Build bamlist
    bamlist: StringVector = []

    # Check if input is from pipe
    if not options.bam_files and sys.stdin.isatty() and not options.bamlist:
        print(TELSEQ_USAGE_MESSAGE)
        sys.exit(1)

    # Process bamlist file if specified
    if options.bamlist:
        bamlist = process_bamlist_from_file(options.bamlist)
    elif not options.bam_files:
        # Read from stdin
        bamlist = process_bamlist_from_stdin()
    else:
        # Use positional arguments
        bamlist = options.bam_files

    if not bamlist:
        print(TELSEQ_USAGE_MESSAGE)
        sys.exit(1)

    # Load exome BED file if specified
    exomebed: Dict[str, List[Range]] = {}
    if options.exomebed:
        exomebed = read_bed_as_dict(options.exomebed)
        print(f"loaded {len(exomebed)} exome regions\n", file=sys.stderr)

    # Create scanner
    scanner = Scanner(
        bam_files=bamlist,
        exome_bed=exomebed if exomebed else None,
        tel_k=options.k,
        ignorerg=options.u,
        mergerg=options.m,
        onebam=options.w,
    )

    # Scan BAM files
    resultlist = scanner.scan()

    # Output results
    output_results(
        resultlist,
        output_file=options.output,
        writerheader=not options.H,
        mergerg=options.m,
    )

    # Print log for exome
    if exomebed:
        printlog(resultlist)

    print("Completed writing results", file=sys.stderr)


def main() -> None:
    """Main entry point."""
    scan_bam()


if __name__ == "__main__":
    main()
