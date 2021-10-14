import glob
import argparse
import pandas as pd
import re
import gff
from os import path
from functools import reduce

parser = argparse.ArgumentParser(description='Tool that combines ijump reports from several files into one summary table')
parser.add_argument('-d', '--dir_report', help='Directory with report files')
parser.add_argument('-o', '--output', help='Output table file')
parser.add_argument('-g', '--gff', nargs='?', default='-', help='Path to gff file')
parser.add_argument('--lab_format',  action='store_true', help='If set, output internal laboratory')
parser.add_argument('--clonal',  action='store_true', help='If set, runs clonal merging')
parser.add_argument('-a', '--a_samples', nargs='?', default='-', help='Path to folder with unevolved population samples. Used for clonal analysis.')
args = parser.parse_args()

report_folder = args.dir_report
out_file = args.output
a_sample_files = []

if args.clonal:
    if args.a_samples == '-':
        print('No unevolved samples are provided')
    else:
        a_sample_files = glob.glob(path.join(args.a_samples, '*.txt'))

report_files = glob.glob(path.join(report_folder, 'ijump_*'))    # collect report files

# Test if there are iJump report files in the directory
if not report_files:
    print('No iJump report files were found. Report files are required to have "ijump_" prefix.')
    exit(1)

# list of samples from reports
sample_list = [re.search('ijump_(.+)\.txt', report).group(1) for report in
               [path.basename(file) for file in report_files]
               ]

# list of unevolved samples
if args.clonal:
    a_sample_list = [re.search('(\dA)\.txt', report).group(1) for report in a_sample_files]

summary_table = pd.DataFrame()

# read all variants
report_dfs = []

for report in report_files:
    sample_name = re.search('ijump_(.+)\.txt', path.basename(report)).group(1)
    print(f'Reading {report}')
    report_dfs.append(pd.read_csv(report, sep='\t').query('Depth > 10').drop(columns='Depth').rename(columns = {'Frequency': sample_name}))

# read unevolved ssmples if provided
if args.clonal:
    for a_sample in a_sample_files:
        sample_name = re.search('(\dA)\.txt', path.basename(a_sample)).group(1)
        print(f'Reading {a_sample}')
        report_dfs.append(pd.read_csv(a_sample, sep='\t').query('Depth > 10').drop(columns='Depth').rename(columns = {'Frequency': sample_name}))

# Merge all data frames to one comparative table
summary_table = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer', on=['IS Name', 'Annotation', 'Chromosome', 'Start', 'Stop']), report_dfs)
summary_table = summary_table.fillna(0)

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
    
    # Add maximum observed frequency
    summary_table['MAX'] = summary_table[sample_list].apply(max, axis=1)
    summary_table = summary_table.query('MAX > 0')
    # List for column names
    cols = []
    cols = ['ID', 'Gene', 'Mutation']
    # Generate short format names
    sample_list_rename = [re.search('_([^_]+)', sample).group(1) for sample in sample_list]

    if args.clonal:
        cols.extend(sample_list_rename)
    else:
        cols.extend(sample_list)

    cols.extend(['MAX', 'Annotation', 'Locus Tag', 'Chromosome'])

    if args.lab_format:
        cols.extend(['MAX', 'Annotation', 'Category', 'Effect', 'Locus Tag', 'Chromosome'])

        if args.clonal:
            # Calculate maximum in unevolved population samples if they are provided
            if args.a_samples != '-':
                summary_table['A_MAX'] = summary_table[a_sample_list].apply(max, axis=1)
                summary_table['Category'] = ['P' if max_a > 0 else 'A' for max_a in
                                                       summary_table['A_MAX'].to_list()]
            else:
                summary_table['A_MAX'] = '?'
                summary_table['Category'] = '?'
            summary_table['Effect'] = 'IS_insertion'
            summary_table['Mutation'] = 'IS_insertion'
            
        else:
            a_samples = [a_sample for a_sample in sample_list if a_sample[-1] == 'A']
            summary_table['A_MAX'] = summary_table[a_samples].apply(max, axis=1)
            summary_table['Effect'] = 'IS_insertion'
            summary_table['Mutation'] = 'IS_insertion'
            summary_table['Category'] = ['P' if max_a > 0 else 'A' for max_a in
                                         summary_table['A_MAX'].to_list()]

        # Rename samples to short format
        summary_table = summary_table.rename(
            columns=dict(zip(sample_list, sample_list_rename))
        )

        if args.clonal:
            if args.a_samples != '-':
                a_sample_list_rename = [re.search('_?([^_]+)', sample).group(1) for sample in a_sample_list]
                summary_table = summary_table.rename(
                    columns=dict(zip(a_sample_list, a_sample_list_rename))
                )
                # Add columns for A-samples to the final table
                cols.extend(a_sample_list)

        cols.extend(['A_MAX'])

        # Sum IS elements with the same name
        summary_table_collapsed = summary_table.drop(columns=['Category']).copy()
        summary_table_collapsed['Mutation'] = [re.match(r'^(.+?)_?\d*$', is_ele).group(1) for
                                               is_ele in summary_table_collapsed['IS Name'].tolist()]
        summary_table_collapsed = summary_table_collapsed.drop(columns=['IS Name', 'Start', 'Stop'])
        summary_table_collapsed = summary_table_collapsed.groupby(
            ['ID', 'Gene', 'Mutation', 'Annotation', 'Effect', 'Locus Tag', 'Chromosome'],
            as_index=False
        ).agg('sum')

        if args.clonal:
            if args.a_samples != '-':
                summary_table_collapsed['Category'] = ['P' if max_a > 0 else 'A' for max_a in summary_table_collapsed['A_MAX'].to_list()]
            else:
                summary_table_collapsed['Category'] = '?'
        else:
            summary_table_collapsed['Category'] = ['P' if max_a > 0 else 'A' for max_a in
                                                   summary_table_collapsed['A_MAX'].to_list()]

        summary_table_collapsed = summary_table_collapsed[cols]

        summary_table_collapsed.to_csv(
            path.join(path.dirname(out_file), 'collapsed_' + path.basename(out_file)),
            sep='\t',
            index=False)

        summary_table_collapsed.query('MAX >= 0.01').to_csv(
            path.join(path.dirname(out_file), 'selected_collapsed_' + path.basename(out_file)),
            sep='\t',
            index=False)

summary_table[cols].to_csv(out_file, sep='\t', index=False)
