#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import re
import subprocess
import sys

import pysam

from ijump import circos
from ijump.isclipped import EstimationMode, ISClipped

# Define a path to output directory that will be available to all functions.
# Required to be global as it will be used to generate output in case of preliminary exit
output_dir = "."


# Build the makeblastdb command as an argv list, kept separate from execution
# so the construction can be unit-tested without invoking BLAST+.
def makeblastdb_command(ref_name, ref_file):
    return ["makeblastdb", "-in", ref_file, "-dbtype", "nucl", "-out", ref_name]


# Check if BLAST database file exists for reference genome. If not -create it.
def check_blast_ref(ref_name, ref_file):
    if os.path.isfile(ref_name + ".nsq"):  # if blast database exists pass or make it for reference
        pass
    else:
        try:
            subprocess.run(
                makeblastdb_command(ref_name, ref_file),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except FileNotFoundError:
            logging.error("makeblastdb not found. Is BLAST+ installed and on PATH?")
            raise
        except subprocess.CalledProcessError as e:
            logging.error(
                f"makeblastdb failed (exit {e.returncode}): {e.stderr.decode(errors='replace')}"
            )
            raise
        if not os.path.isfile(ref_name + ".nsq"):
            raise RuntimeError(f"makeblastdb reported success but {ref_name}.nsq was not created")


def build_arg_parser():
    # Command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "iJump searches for small frequency IS elements rearrangements in evolved populations"
        )
    )
    parser.add_argument("-a", "--aln", type=str, action="store", help="BAM or SAM alignment file")
    parser.add_argument(
        "-r", "--ref", type=str, action="store", help="Reference genome in FASTA format"
    )
    parser.add_argument(
        "-g",
        "--gff",
        type=str,
        action="store",
        help="Annotations in GFF format for reference genome. Required for average mode.",
    )
    parser.add_argument(
        "-i", "--isel", type=str, action="store", help="File with IS elements coordinates"
    )
    parser.add_argument(
        "-c",
        "--circos",
        action="store_true",
        default=False,
        help="Set flag to build input files for CIRCOS",
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=".", help="Output directory. Default: . (current)"
    )
    parser.add_argument(
        "-w",
        "--wd",
        type=str,
        default="ijump_wd",
        help="Work directory. Default: ijump_wd (current)",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=200,
        help="Radius around IS elements boundaries to search soft clipped reads.",
    )
    # default=EstimationMode.AVERAGE.value (the plain string "average"), not
    # EstimationMode.AVERAGE itself: argparse runs `type` over the default
    # whenever it is a str instance, and EstimationMode members are str
    # instances too (str, Enum), so a default of EstimationMode.AVERAGE would
    # go through type=str and become str(EstimationMode.AVERAGE) ==
    # "EstimationMode.AVERAGE" -- a value that matches neither
    # EstimationMode.AVERAGE nor EstimationMode.PRECISE and silently no-ops
    # ISClipped.run()'s mode dispatch. Using the plain value string sidesteps
    # that: type=str leaves it unchanged. parse_args() below converts the
    # resulting plain string (default or explicitly typed) into the real
    # EstimationMode member, uniformly for both paths.
    parser.add_argument(
        "--estimation_mode",
        type=str,
        default=EstimationMode.AVERAGE.value,
        choices=list(EstimationMode),
        help="Specifies how the IS frequency will be esimated. 'average' - by averaging"
        " the region coverage and number of clipped reads. Or 'precise' - iJump will"
        " try to separate each insertion event.",
    )
    parser.add_argument("--version", action="store_true", help="Print iJump version and exit.")
    return parser


def parse_args(argv=None):
    args = build_arg_parser().parse_args(argv)
    # Normalize to a genuine EstimationMode member: argparse's choices check
    # above ran against the plain string (so an invalid value still reports
    # "invalid choice"), and args.estimation_mode is a plain str at this
    # point regardless of whether it came from the default or an explicit
    # --estimation_mode flag.
    args.estimation_mode = EstimationMode(args.estimation_mode)
    return args


def main():
    args = parse_args()

    global output_dir
    output_dir = args.outdir

    # Make output directory if not exists.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize logger.
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    log_formatter = logging.Formatter("%(levelname)s: %(asctime)s - %(message)s")

    file_handler = logging.FileHandler(os.path.join(output_dir, "ijump.log"))
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    # Print iJump version.
    version = "1.0.4"
    logging.info(f"iJump v.{version}\n")
    logging.info("author: Semion Leyn")
    logging.info("Please ask questions and report issues on GitHub page of the project:")
    logging.info("https://github.com/sleyn/ijump\n")
    logging.info(f"Mode: {args.estimation_mode}")
    logging.info(f"Alignment file: {args.aln}\n")

    # Collect required file information from arguments.
    alignment_file = args.aln
    gff = args.gff
    is_file = args.isel
    reference_filename = re.match(r"(.+)\.", args.ref)
    reference = reference_filename.group(1)

    # Check exisitance of  BLAST database.
    check_blast_ref(reference, args.ref)

    # Check alignment file type (SAM/BAM).
    if alignment_file[-3:] == "sam":
        a_type = ""
    elif alignment_file[-3:] == "bam":
        a_type = "b"

    # Make work directory if not exists
    if not os.path.exists(args.wd):
        os.makedirs(args.wd)
    else:
        if len(os.listdir(args.wd)) != 0:
            logging.warning("The work directory is not empty. Will be cleanned.")
            # Clean work directory
            for file in glob.glob(os.path.join(args.wd, "*")):
                os.remove(file)

    # Open read alignment file (SAM or BAM)
    alignment = pysam.AlignmentFile(alignment_file, "r" + a_type)

    # Process alignment to estimate IS elements frequencies.
    # Initialize new instance of ISClipped class object.
    is_processing = ISClipped(
        alignment, reference, gff, args.wd, os.path.join(output_dir, "ijump_junction_pairs.txt")
    )

    # Collect IS coordinates from IS file.
    is_processing.iscollect(is_file)

    # Set area to search for clipped reads.
    is_processing.set_is_boundaries(args.radius)

    # Run the average/precise pipeline and write its output files.
    result = is_processing.run(args.estimation_mode)
    if not result.insertions_found:
        logging.info(result.message)
        sys.exit(0)

    # Plot circular diagram of insertions
    if args.circos is True and args.estimation_mode == EstimationMode.AVERAGE:
        circos.write_files(
            is_processing.report_table,
            is_processing.sum_by_region,
            is_processing.is_coords,
            is_processing.ref_len,
            is_processing.data_folder,
            is_processing.cutoff,
            is_processing.average_depth,
            is_processing.gff.ann_pos,
        )


if __name__ == "__main__":
    main()
