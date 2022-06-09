import os
import argparse

import pandas as pd
import pysam
import isclipped
import re
import subprocess

# Check if junctions exist for IS elements
def check_junctions_presence(junc_tbl, outdir):
    if junc_tbl.size:
        junc_tbl.to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep='\t', index=False)
    else:
        print('No junctions was found')
        exit(0)


if __name__ == "__main__":
    # read arguments
    parser = argparse.ArgumentParser(description="iJump searches for small frequency IS elements rearrangements in evolved populations")
    parser.add_argument('-a', '--aln', type=str, action='store', help='BAM or SAM alignment')
    parser.add_argument('-r', '--ref', type=str, action='store', help='Reference genome in FASTA format')
    parser.add_argument('-g', '--gff', type=str, action='store', help='Annotations in GFF format for reference genome')
    parser.add_argument('-i', '--isel', type=str, action='store', help='File with IS elements coordinates')
    parser.add_argument('-c', '--circos', action='store_true', default=False, help="Set flag to build input files for CIRCOS")
    parser.add_argument('-o', '--outdir', type=str, default='.', help="Output directory")
    parser.add_argument('--radius', type=int, default=200, help="Radius around IS elements boundaries to search soft clipped reads.")
    parser.add_argument('--estimation_mode', type=str, default='average',
                        help="Specifies how the IS frequency will be esimated. 'average' - by averaging the region coverage"
                             " and number of clipped reads. Or 'point' - iJump will try to separate each insertion event.")
    args = parser.parse_args()

    alignment_file = args.aln
    gff = args.gff
    is_file = args.isel
    filename = re.match('(.+)\.', args.ref)
    reference = filename.group(1)

    if os.path.isfile(reference + '.nsq'):          # if blast database exists pass or make it for reference
        pass
    else:
        makeblastdb_command = "makeblastdb -in " + args.ref + " -dbtype nucl -out " + reference
        makeblastdb = subprocess.Popen(makeblastdb_command.split(), stdout=subprocess.PIPE)

    # check alignment type
    if alignment_file[-3:] == 'sam':
        a_type = ''
    elif alignment_file[-3:] == 'bam':
        a_type = 'b'

    # Make output directories if not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)


    # read alignment
    alignment = pysam.AlignmentFile(alignment_file, "r" + a_type)

    x = isclipped.isclipped(alignment, reference, gff)
    x.outdir = args.outdir
    x.gff.readgff()
    x.gff.pos_to_ann()
    x.iscollect(is_file)
    x.crtable(args.radius)   # provide radius
    x.clipped_reads.to_csv(os.path.join(args.outdir, "reads.txt"), sep='\t', index=False)
    x.runblast('cl.fasta', 'cl_blast.out', 1)
    x.parseblast(1, 'cl_blast.out')
    if args.estimation_mode == 'average':
        x.call_junctions(1)
        check_junctions_presence(x.junctions, args.outdir)
        x.summary_junctions_by_region()
        x.sum_by_region.to_csv(os.path.join(args.outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)
        x.report()
        x.report_table.to_csv(os.path.join(args.outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)
    elif args.estimation_mode == 'point':
        # Make table of regions in the reference genome where extract clipped reads for backwards assignment
        reference_regions = x.make_gene_side_regions()
        # We will need average read length to extend region boudaries
        x.av_read_len = x.read_lengths / x.n_reads_analyzed
        x.crtable_bwds(reference_regions)
        x.runblast('cl_bwrd.fasta', 'cl_blast_bwds.out', 0)
        x.parseblast(0, 'cl_blast_bwds.out')
        x.call_junctions(0)
        check_junctions_presence(x.junctions, args.outdir)
        x.search_insert_pos()

        positions = pd.concat(
            [
                x.pairs_df[['Position_l', 'Chrom']].rename(columns={'Position_l': 'Position'}),
                x.pairs_df[['Position_r', 'Chrom']].rename(columns={'Position_r': 'Position'})
            ],
            axis=0
        ).drop_duplicates()
        x.count_depth(positions)
        x.assess_isel_freq()
        x.pairs_df.to_csv(os.path.join(args.outdir, "ijump_junction_pairs.txt"), sep='\t', index=False)

    if args.circos is True:
        x.create_circos_files()