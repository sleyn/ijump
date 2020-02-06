import glob
import argparse
import pandas as pd
import sys
import re
import gff

parser = argparse.ArgumentParser(description='Tool that combines ijump reports from several files into one summary table')
parser.add_argument('-d', '--dir_report', help='Directory with report files')
parser.add_argument('-o', '--output', help='Output table file')
parser.add_argument('-g', '--gff', nargs='?', default='-', help='Path to gff file')
args = parser.parse_args()

report_folder = args.dir_report
out_file = args.output

report_files = glob.glob(report_folder + '/ijump_*')    # collect report files

# Test if there are iJump report files in the directory
if not report_files:
    print('No iJump report files were found. Report files are required to have "ijump_" prefix.')
    exit(1)

sample_list = list()            # list of samples from reports
summary_table = pd.DataFrame()

# collect all variants
for report in report_files:
    print(f'Reading {report}')
    match = re.search('[-|_](.+)\.txt', report)
    sample_list.append(match.group(1))

    report_df = pd.read_csv(report,  header=0, sep='\t')
    summary_table = summary_table.append(report_df[['Chromosome', 'Annotation', 'Start', 'Stop', 'IS Name']], sort=False)

summary_table = summary_table.drop_duplicates()         # remove duplicated strings
st_len = len(summary_table)                             # length of summary table
sample_list.sort()
for sample in sample_list:                              # fill columns of summary table with 0s
    summary_table[sample] = [0 for x in range(st_len)]

summary_table['index'] = [str(summary_table.iloc[i]['Start']) + summary_table.iloc[i]['IS Name'] + str(summary_table.iloc[i]['Stop']) for i in range(len(summary_table))]
summary_table.set_index('index', inplace=True)
summary_table[sample_list] = summary_table[sample_list].astype('float')

for report in report_files:
    report_df = pd.read_csv(report, header=0, sep='\t')
    match = re.search('[-|_](.+)\.txt', report)
    sample_name = match.group(1)

    for i in range(len(report_df)):
        if report_df.iloc[i]['Depth'] > 10:
            ind = str(report_df.iloc[i]['Start']) + report_df.iloc[i]['IS Name'] + str(report_df.iloc[i]['Stop'])
            summary_table.at[ind, sample_name] = report_df.iloc[i]['Frequency']

if args.gff == '-':
    summary_table['Functional annotation'] = ['-'] * len(summary_table)
else:
    gff_file = gff.gff(args.gff)
    gff_file.readgff()
    summary_table['Functional annotation'] = summary_table.apply(lambda x: gff_file.gff_pos[x['Chromosome']][int( (int(x['Start']) + int(x['Stop']) ) / 2)][3], axis=1)

summary_table.to_csv(out_file, sep='\t', index=False)

