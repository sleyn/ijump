import glob
import argparse
import pandas as pd
import re
import gff
from os import path
from functools import reduce


# Read and ijump report tables
def collect_reports(report_path, a_samples_path, clonal_workflow=False):
    # Collect junctions from unevolved samples
    a_sample_files = []

    if clonal_workflow:
        if a_samples_path == '-':
            print('No unevolved samples are provided')
        else:
            a_sample_files = glob.glob(path.join(a_samples_path, '*.txt'))

    # Collect paths to report tables for samples.
    report_files = glob.glob(path.join(report_path, 'ijump_*'))

    # Test if there are iJump report files in the directory
    if not report_files:
        print('No iJump report files were found. Report files are required to have "ijump_" prefix.')
        exit(1)

    # List of samples from reports
    sample_list = [re.search('ijump_(.+)\.txt', report).group(1) for report in
                   [path.basename(file) for file in report_files]
                   ]
    sample_list.sort()

    a_sample_list = []

    # List of unevolved samples
    if clonal_workflow:
        a_sample_list = [re.search('(\dA)\.txt', report).group(1) for report in a_sample_files]

    return report_files, a_sample_files, sample_list, a_sample_list


def read_reports(report_files, a_sample_files, clonal_workflow, mode):
    # Read all variants
    report_dfs = []

    for report in report_files:
        sample_name = re.search('ijump_(.+)\.txt', path.basename(report)).group(1)
        print(f'Reading {report}')
        if mode == 'average':
            report_dfs.append(
                pd.read_csv(report, sep='\t').
                    query('Depth > 10').
                    drop(columns='Depth').
                    rename(columns={'Frequency': sample_name})
            )
        elif mode == 'precise':
            report_dfs.append(
                pd.read_csv(report, sep='\t')[
                    ['Position_l', 'Position_r', 'Chrom', 'IS_name', 'Depth', 'Frequency']].
                    query('Depth > 10').
                    drop(columns='Depth').
                    rename(
                        columns={
                            'Frequency': sample_name,
                            'Position_l': 'Start',
                            'Position_r': 'Stop',
                            'IS_name': 'IS Name',
                            'Chrom': 'Chromosome'
                        }
                    )
            )

    # Read unevolved samples if provided (only for clonal data)
    if clonal_workflow:
        for a_sample in a_sample_files:
            sample_name = re.search(r'(\dA)\.txt', path.basename(a_sample)).group(1)
            print(f'Reading {a_sample}')
            if mode == 'average':
                report_dfs.append(pd.read_csv(a_sample, sep='\t').query('Depth > 10').drop(columns='Depth').rename(
                    columns={'Frequency': sample_name}))
            elif mode == 'precise':
                report_dfs.append(pd.read_csv(a_sample, sep='\t')[
                    ['Position_l', 'Position_r', 'Chrom', 'IS_name', 'Depth', 'Frequency']].
                    query('Depth > 10').
                    drop(columns='Depth').
                    rename(
                        columns={
                            'Frequency': sample_name,
                            'Position_l': 'Start',
                            'Position_r': 'Stop',
                            'IS_name': 'IS Name',
                            'Chrom': 'Chromosome'
                        }
                    )
                )

    # Merge all data frames to one comparative table
    if mode == 'average':
        id_columns = ['IS Name', 'Annotation', 'Chromosome', 'Start', 'Stop']
    elif mode == 'precise':
        id_columns = ['IS Name', 'Chromosome', 'Start', 'Stop']

    summary_table = reduce(
        lambda df1, df2: pd.merge(df1, df2, how='outer', on=id_columns),
        report_dfs)
    summary_table = summary_table.fillna(0)
    return summary_table


# Annotate junctions
def annotate_feild(gff_file, chrom, pos_1, pos_2, field, mode):
    annotation_fields = {
        'Locus Tag': 0,
        'Gene': 1,
        'ID': 2,
        'Annotation': 3
    }

    field_ix = annotation_fields[field]

    if mode == 'average':
        return gff_file.gff_pos[chrom][int((int(pos_1) + int(pos_2)) / 2)][field_ix]
    elif mode == 'precise':
        l_ann = gff_file.gff_pos[chrom][int(pos_1)][field_ix]
        r_ann = gff_file.gff_pos[chrom][int(pos_1)][field_ix]
        if l_ann == r_ann:
            return l_ann
        else:
            return f'{l_ann} <> {r_ann}'


# Add GFF annotations.
def add_gff_annotations(summary_table, gff_file_path, mode):
    # If no gff file was provided just put '-' in all annotation fields.
    if gff_file_path == '-':
        n_variants = len(summary_table)
        summary_table['Locus Tag'] = ['-'] * n_variants
        summary_table['Gene'] = ['-'] * n_variants
        summary_table['ID'] = ['-'] * n_variants
        summary_table['Annotation'] = ['-'] * n_variants
    else:
        # Read GFF file.
        gff_file = gff.gff(gff_file_path)
        gff_file.readgff()

        # Add annotations.
        for field in ['Locus Tag', 'Gene', 'ID', 'Annotation']:
            summary_table[field] = summary_table.apply(
                lambda x: annotate_feild(gff_file, x['Chromosome'], x['Start'], x['Stop'], field, mode),
                axis=1)

    return summary_table


def main():
    parser = argparse.ArgumentParser(
        description='Tool that combines ijump reports from several files into one summary table')
    parser.add_argument('-d', '--dir_report', help='Directory with report files')
    parser.add_argument('-o', '--output', help='Output table file')
    parser.add_argument('-g', '--gff', nargs='?', default='-', help='Path to gff file')
    parser.add_argument('-p', '--prefix', default='', help='If set would be used as prefix')
    parser.add_argument('-m', '--ijump_mode', default='average', help='Variant of iiJump pipeline: "average" or'
                                                                      ' "precise". Default: "average"')
    parser.add_argument('--lab_format', action='store_true', help='If set, output internal laboratory')
    parser.add_argument('--clonal', action='store_true', help='If set, runs clonal merging')
    parser.add_argument('-a', '--a_samples', nargs='?', default='-',
                        help='Path to folder with unevolved population samples. Used for clonal analysis.')
    args = parser.parse_args()

    # Check if correct mode was given.
    if args.ijump_mode not in ['average', 'precise']:
        print('Wrong iJump mode. Choose "average" or "precise".')
        exit(0)

    # Folder with ijump report files
    report_dir = args.dir_report
    # Path to output table
    out_file = args.output

    # Collect report files for evolved and unevolved samples.
    report_files, a_sample_files, sample_list, a_sample_list = collect_reports(report_dir, args.a_samples, args.clonal)
    # Read all reports and combine them to the summary table
    summary_table = read_reports(report_files, a_sample_files, args.clonal, args.ijump_mode)
    # Add annotations from GFF file
    summary_table = add_gff_annotations(summary_table, args.gff, args.ijump_mode)

    if args.gff != '-':
        # Add maximum observed frequency.
        summary_table['MAX'] = summary_table[sample_list].apply(max, axis=1)
        summary_table = summary_table.query('MAX > 0').copy()

        # Initiate list for column names.
        cols = ['Start', 'Stop', 'ID', 'Gene', 'Mutation']

        if not args.lab_format:
            cols.extend(sample_list)
            cols.extend(['MAX', 'Annotation', 'Locus Tag', 'Chromosome'])
        else:
            # Generate short format names.
            sample_list_rename = [re.search('_([^_]+)$', sample).group(1) for sample in sample_list]
            cols.extend(sample_list_rename)
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
                summary_table['Mutation'] = summary_table['IS Name']
            else:
                # Process population unevolved samples
                a_pop_samples = [a_sample for a_sample in sample_list if a_sample[-1] == 'A']
                # If unevolved population samples are provided calculate max frequency in unevolved
                if len(a_pop_samples) > 0:
                    summary_table['A_MAX'] = summary_table[a_pop_samples].apply(max, axis=1)
                    summary_table['Category'] = ['P' if max_a > 0 else 'A' for max_a in
                                                 summary_table['A_MAX'].to_list()]
                else:
                    summary_table['A_MAX'] = '?'
                    summary_table['Category'] = '?'

                summary_table['Effect'] = 'IS_insertion'
                summary_table['Mutation'] = summary_table['IS Name']

            # Rename samples to short format
            summary_table = summary_table.rename(
                columns=dict(zip(sample_list, sample_list_rename))
            )

            if args.clonal and args.a_samples != '-':
                a_sample_list_rename = [re.search('_?([^_]+)$', sample).group(1) for sample in a_sample_list]
                summary_table = summary_table.rename(
                    columns=dict(zip(a_sample_list, a_sample_list_rename))
                )
                # Add columns for A-samples to the final table
                cols.extend(a_sample_list)

            cols.extend(['A_MAX'])

            # Sum IS elements with the same name
            summary_table_collapsed = summary_table.copy()
            if args.ijump_mode == 'average':
                summary_table_collapsed['Mutation'] = [re.match(r'^(.+?)_?\d*$', is_ele).group(1) for
                                                       is_ele in summary_table_collapsed['IS Name'].tolist()]
                summary_table_collapsed = summary_table_collapsed.drop(columns=['IS Name', 'Start', 'Stop'])
                summary_table_collapsed = summary_table_collapsed.groupby(
                    ['ID', 'Gene', 'Mutation', 'Annotation', 'Effect', 'Locus Tag', 'Chromosome'],
                    as_index=False
                ).agg('sum')

                if args.a_samples != '-' or len(a_pop_samples) > 0:
                    summary_table_collapsed['Category'] = ['P' if max_a > 0 else 'A' for max_a in
                                                           summary_table_collapsed['A_MAX'].to_list()]
                else:
                    summary_table_collapsed['Category'] = '?'
                    summary_table_collapsed['A_MAX'] = '?'

                summary_table_collapsed = summary_table_collapsed[cols]

                # Write collapsed table. All ISname_"[digit]" frequences will be sum to one ISname field frequency.
                out_file_collapsed = '_'.join(
                    [part for part in [args.prefix, 'collapsed', path.basename(out_file)] if part])
                summary_table_collapsed.to_csv(
                    path.join(path.dirname(out_file), out_file_collapsed),
                    sep='\t',
                    index=False)

            # Write selected (Frequency >= 1%) IS insertions variant table to file.
            out_file_sel_collapsed = '_'.join(
                [part for part in [args.prefix, 'selected_collapsed', path.basename(out_file)] if part])
            summary_table_collapsed.query('MAX >= 0.01').to_csv(
                path.join(path.dirname(out_file), out_file_sel_collapsed),
                sep='\t',
                index=False)

        # Write full IS insertions variant table to file.
        out_file_full = '_'.join([part for part in [args.prefix, path.basename(out_file)] if part])
        summary_table[cols].to_csv(
            path.join(
                path.dirname(out_file),
                out_file_full),
            sep='\t',
            index=False
        )


if __name__ == "__main__":
    main()
