import os
import glob
import argparse
import pandas as pd
import pysam
from isclipped import ISClipped
import re
import subprocess
import logging
import sys


# Check if junctions exist for IS elements
def check_junctions_presence(junc_tbl, outdir):
    if junc_tbl.size:
        junc_tbl.to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep='\t', index=False)
    else:
        logging.info('No junctions was found')
        exit(0)


# Check if BLAST database file exists for reference genome. If not -create it.
def check_blast_ref(ref_name, ref_file):
    if os.path.isfile(ref_name + '.nsq'):  # if blast database exists pass or make it for reference
        pass
    else:
        makeblastdb_command = f'makeblastdb -in {ref_file} -dbtype nucl -out {ref_name}'
        makeblastdb = subprocess.Popen(makeblastdb_command.split(), stdout=subprocess.PIPE)


def main():
    # Command line arguments
    parser = argparse.ArgumentParser(
        description="iJump searches for small frequency IS elements rearrangements in evolved populations")
    parser.add_argument('-a', '--aln', type=str, action='store', help='BAM or SAM alignment file')
    parser.add_argument('-r', '--ref', type=str, action='store', help='Reference genome in FASTA format')
    parser.add_argument('-g', '--gff', type=str, action='store', help='Annotations in GFF format for reference genome. '
                                                                      'Required for average mode.')
    parser.add_argument('-i', '--isel', type=str, action='store', help='File with IS elements coordinates')
    parser.add_argument('-c', '--circos', action='store_true', default=False,
                        help="Set flag to build input files for CIRCOS")
    parser.add_argument('-o', '--outdir', type=str, default='.', help="Output directory. Default: . (current)")
    parser.add_argument('-w', '--wd', type=str, default='ijump_wd', help="Work directory. Default: ijump_wd (current)")
    parser.add_argument('--radius', type=int, default=200,
                        help="Radius around IS elements boundaries to search soft clipped reads.")
    parser.add_argument('--estimation_mode', type=str, default='average',
                        help="Specifies how the IS frequency will be esimated. 'average' - by averaging the region coverage"
                             " and number of clipped reads. Or 'precise' - iJump will try to separate each insertion event.")
    parser.add_argument('--version', action='store_true', help='Print iJump version and exit.')
    args = parser.parse_args()

    # Initialize logger.
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    log_formatter = logging.Formatter("%(levelname)s: %(asctime)s - %(message)s")

    file_handler = logging.FileHandler(os.path.join(args.outdir, 'ijump.log'))
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    # Print iJump version.
    version = '1.0'
    logging.info(f'iJump v.{version}\n')
    logging.info(f'author: Semion Leyn')
    logging.info(f'Please ask questions and report issues on GitHub page of the project:')
    logging.info(f'https://github.com/sleyn/ijump\n')
    logging.info(f'Mode: {args.estimation_mode}')
    logging.info(f'Alignment file: {args.aln}\n')

    # Collect required file information from arguments.
    alignment_file = args.aln
    gff = args.gff
    is_file = args.isel
    reference_filename = re.match(r'(.+)\.', args.ref)
    reference = reference_filename.group(1)

    # Check exisitance of  BLAST database.
    check_blast_ref(reference, args.ref)

    # Check alignment file type (SAM/BAM).
    if alignment_file[-3:] == 'sam':
        a_type = ''
    elif alignment_file[-3:] == 'bam':
        a_type = 'b'

    # Make work directory if not exists
    if not os.path.exists(args.wd):
        os.makedirs(args.wd)
    else:
        if len(os.listdir(args.wd)) != 0:
            logging.warning('The work directory is not empty. Will be cleanned.')
            # Clean work directory
            for file in glob.glob(os.path.join(args.wd, '*')):
                os.remove(file)

    # Make output directory if not exists.
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Open read alignment file (SAM or BAM)
    alignment = pysam.AlignmentFile(alignment_file, "r" + a_type)

    # Process alignment to estimate IS elements frequencies.
    # Initialize new instance of ISClipped class object.
    is_processing = ISClipped(alignment, reference, gff, args.wd)

    # Collect IS coordinates from IS file.
    is_processing.iscollect(is_file)

    # Set area to search for clipped reads.
    is_processing.set_is_boundaries(args.radius)

    # Collect clipped reads inforamtion.
    is_processing.collect_clipped_reads()
    is_processing.clipped_reads.to_csv(os.path.join(args.outdir, "reads.txt"), sep='\t', index=False)

    # Run BLAST to search insertion positions in Reference.
    # 1 - search in IS->Reference direction.
    is_processing.runblast('cl.fasta', 'cl_blast.out', 1)
    is_processing.parseblast(1, 'cl_blast.out')

    # Read GFF file.
    is_processing.gff.readgff()
    is_processing.gff.pos_to_ann()

    # Workflow in "average" mode
    if args.estimation_mode == 'average':
        # Make a table of observed junction positions
        is_processing.call_junctions(1)

        # Check if any junction is present. If not - stop the workflow.
        check_junctions_presence(is_processing.junctions, args.outdir)

        # Count reads supporting IS elements insertions for each IS element and each GE
        is_processing.sum_by_reg_tbl_init()
        is_processing.summary_junctions_by_region()
        is_processing.sum_by_region.to_csv(os.path.join(args.outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)

        # Make a report table of assessed insertion frequencies in each GE
        is_processing.report_table_init()
        is_processing.report_average()
        is_processing.report_table.to_csv(os.path.join(args.outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)
    elif args.estimation_mode == 'precise':
        # Make table of regions in the reference genome where extract clipped reads for backwards assignment
        reference_regions = is_processing.make_gene_side_regions()

        # We will need average read length to extend region boundaries
        is_processing.av_read_len = is_processing.read_lengths / is_processing.n_reads_analyzed

        # Collect clipped reads at the insertion positions
        # found during forward (IS->Reference) search.
        is_processing.crtable_bwds(reference_regions)

        # Run BLAST to search for positions of found reads.
        is_processing.runblast('cl_bwrd.fasta', 'cl_blast_bwds.out', 0)

        # Collect BLAST results.
        is_processing.parseblast(0, 'cl_blast_bwds.out')

        # Format results as junction table to attribute reads and their junction positions to the IS elements.
        is_processing.call_junctions(0)

        # Check if any junction is present. If not - stop the workflow.
        check_junctions_presence(is_processing.junctions, args.outdir)

        # Find pairs of junctions that should indicate insertion positions of both edges of IS element.
        is_processing.search_insert_pos()

        positions = pd.concat(
            [
                is_processing.pairs_df[['Position_l', 'Chrom']].rename(columns={'Position_l': 'Position'}),
                is_processing.pairs_df[['Position_r', 'Chrom']].rename(columns={'Position_r': 'Position'})
            ],
            axis=0
        ).drop_duplicates()

        # Count depth of unclipped reads to have a background depth of coverage
        is_processing.count_depth_unclipped(positions)

        # Make an estimate of insertion frequency
        is_processing.assess_isel_freq()
        is_processing.pairs_df.to_csv(os.path.join(args.outdir, "ijump_junction_pairs.txt"), sep='\t', index=False)

    # Plot circular diagram of insertions
    if args.circos is True and args.estimation_mode == 'precise':
        is_processing.create_circos_files()


if __name__ == "__main__":
    main()
