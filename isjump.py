import sys
import os
import getopt
import pysam
import isclipped
import re
import subprocess

# read arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:r:g:i:")
except getopt.GetoptError:
    print("isfinder_parse.py -f <Alignment file> -r <Reference in a fasta format> -g <GFF file> -i <IS_File>")
    sys.exit(1)

for opt, arg in opts:
    if opt == '-f':
        alignment_file = arg
    elif opt == '-r':
        filename = re.match('(.+)\.', arg)
        reference = filename.group(1)
        if os.path.isfile(reference + '.nsq'):          # if blast database exists pass or make it for reference
            pass
        else:
            makeblastdb_command = "makeblastdb -in " + arg + " -dbtype nucl -out " + reference
            makeblastdb = subprocess.Popen(makeblastdb_command.split(), stdout=subprocess.PIPE)
    elif opt == '-g':
        gff = arg
    elif opt == '-i':
        is_file = arg
    else:
        print("isfinder_parse.py -f <Alignment file>")
        sys.exit(1)

# check alignment type
if alignment_file[-3:] == 'sam':
    a_type = ''
elif alignment_file[-3:] == 'bam':
    a_type = 'b'

# read IS coordinares from file


# read alignment
alignment = pysam.AlignmentFile(alignment_file, "r" + a_type)

x = isclipped.isclipped(alignment, reference, gff)
x.gff.readgff()
x.gff.pos_to_ann()
x.iscollect(is_file)
x.crtable(200)   # provide radius
x.clipped_reads.to_csv("reads.txt", sep='\t', index=False)
x.runblast()
x.parseblast()
x.call_junctions()
x.junctions.to_csv("isjump_junctions.txt", sep='\t', index=False)
x.summary_junctions_by_region()
x.sum_by_region.to_csv("isjump_sum_by_reg.txt", sep='\t', index=False)
x.report()
x.report_table.to_csv("isjump_report_by_is_reg.txt", sep='\t', index=False)
x.create_circos_files()
