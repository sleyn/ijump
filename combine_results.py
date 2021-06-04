import glob
import argparse
import pandas as pd
import re
import gff
from os import path

parser = argparse.ArgumentParser(description='Tool that combines ijump reports from several files into one summary table')
parser.add_argument('-d', '--dir_report', help='Directory with report files')
parser.add_argument('-o', '--output', help='Output table file')
parser.add_argument('-g', '--gff', nargs='?', default='-', help='Path to gff file')
parser.add_argument('--lab_format',  action='store_true', help='If set, output internal laboratory')
parser.add_argument('--clonal',  action='store_true', help='If set, runs clonal merging')
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
    summary_table = summary_table.append(report_df[['Chromosome', 'Start', 'Stop', 'IS Name']], sort=False)

summary_table = summary_table.drop_duplicates()         # remove duplicated strings
st_len = len(summary_table)                             # length of summary table
sample_list.sort()
for sample in sample_list:                              # fill columns of summary table with 0s
    summary_table[sample] = [0 for x in range(st_len)]

summary_table['index'] = [str(summary_table.iloc[i]['Start']) + summary_table.iloc[i]['IS Name'] + str(summary_table.iloc[i]['Stop']) for i in range(len(summary_table))]
summary_table.set_index('index', inplace=True)
summary_table[sample_list] = summary_table[sample_list].astype('float')

for report in report_files:
    # read report file and extract sample name
    report_df = pd.read_csv(report, header=0, sep='\t')
    match = re.search('[-|_](.+)\.txt', report)
    sample_name = match.group(1)

    for i in range(len(report_df)):
        if report_df.iloc[i]['Depth'] > 10:
            ind = str(report_df.iloc[i]['Start']) + report_df.iloc[i]['IS Name'] + str(report_df.iloc[i]['Stop'])
            summary_table.at[ind, sample_name] = report_df.iloc[i]['Frequency']

# if no gff file was provided just put '-' in all annotation fields
if args.gff == '-':
    summary_table['Locus Tag'] = ['-'] * len(summary_table)
    summary_table['Gene'] = ['-'] * len(summary_table)
    summary_table['ID'] = ['-'] * len(summary_table)
    summary_table['Annotation'] = ['-'] * len(summary_table)
else:
    gff_file = gff.gff(args.gff)
    gff_file.readgff()
    summary_table['Locus Tag'] = summary_table.apply(lambda x: gff_file.gff_pos[x['Chromosome']][int( (int(x['Start']) + int(x['Stop']) ) / 2)][0], axis=1)
    summary_table['Gene'] = summary_table.apply(lambda x: gff_file.gff_pos[x['Chromosome']][int( (int(x['Start']) + int(x['Stop']) ) / 2)][1], axis=1)
    summary_table['ID'] = summary_table.apply(lambda x: gff_file.gff_pos[x['Chromosome']][int( (int(x['Start']) + int(x['Stop']) ) / 2)][2], axis=1)
    summary_table['Annotation'] = summary_table.apply(lambda x: gff_file.gff_pos[x['Chromosome']][int( (int(x['Start']) + int(x['Stop']) ) / 2)][3], axis=1)
    
    # rearrange columns so all annotation columns should go first
    cols = summary_table.columns.tolist()
    cols = cols[:4] + cols[-4:] + cols[4:-4]
    summary_table = summary_table[cols]
    summary_table['MAX'] = summary_table[sample_list].apply(max, axis=1)

    if args.lab_format:
        if not args.clonal:
            a_samples = [a_sample for a_sample in sample_list if a_sample[-1] == 'A']
            summary_table['A_MAX'] = summary_table[a_samples].apply(max, axis=1)
            summary_table['Category'] = ['P' if max_a > 0 else 'A' for max_a in summary_table['A_MAX'].to_list()]
            summary_table['Effect'] = 'IS_insertion'
            summary_table['Mutation'] = 'IS_insertion'
            cols = ['ID', 'Gene', 'Mutation']
            cols.extend(sample_list)
            cols.extend(['MAX', 'A_MAX', 'Annotation', 'Category', 'Effect', 'Locus Tag', 'Chromosome', 'IS Name'])
            summary_table = summary_table[cols]

            # Sum IS elements with the same name
            summary_table_collapsed = summary_table.copy()
            summary_table_collapsed['Mutation'] = [re.match(r'^(.+?)_?\d*$', is_ele).group(1) for
                                                 is_ele in summary_table_collapsed['IS Name'].tolist()]
            summary_table_collapsed = summary_table_collapsed.drop(columns=['IS Name'])
            summary_table_collapsed.groupby(
                ['ID', 'Gene', 'Mutation', 'Annotation', 'Category', 'Effect', 'Locus Tag', 'Chromosome'],
                as_index = False
            ).agg('sum')

            summary_table_collapsed.to_csv(
                path.join(path.dirname(out_file), 'collapsed' + path.basename(out_file)),
                sep='\t',
                index=False)


summary_table.to_csv(out_file, sep='\t', index=False)
