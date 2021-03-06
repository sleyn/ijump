import os
import argparse
import pysam
import isclipped
import re
import subprocess

# read arguments
parser = argparse.ArgumentParser(description="iJump searches for small frequency IS elements rearrangements in evolved populations")
parser.add_argument('-a', '--aln', type=str, action='store', help='BAM or SAM alignment')
parser.add_argument('-r', '--ref', type=str, action='store', help='Reference genome in FASTA format')
parser.add_argument('-g', '--gff', type=str, action='store', help='Annotations in GFF format for reference genome')
parser.add_argument('-i', '--isel', type=str, action='store', help='File with IS elements coordinates')
parser.add_argument('-c', '--circos', action='store_true', default=False, help="Set flag to build input files for CIRCOS")
parser.add_argument('-o', '--outdir', type=str, default='.', help="Output directory")
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
x.crtable(200)   # provide radius
x.clipped_reads.to_csv(os.path.join(args.outdir, "reads.txt"), sep='\t', index=False)
x.runblast()
x.parseblast()
x.call_junctions()
if x.junctions.size:
    x.junctions.to_csv(os.path.join(args.outdir, "ijump_junctions.txt"), sep='\t', index=False)
else:
    print('No junctions was found')
    exit(0)
x.summary_junctions_by_region()
x.sum_by_region.to_csv(os.path.join(args.outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)
x.report()
x.report_table.to_csv(os.path.join(args.outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)

if args.circos is True:
    x.create_circos_files()