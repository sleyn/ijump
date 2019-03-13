import glob
import getopt
import pandas as pd
import sys
import re

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:")
except getopt.GetoptError:
    print("combine_results.py -d <Directory with report files>")
    sys.exit(1)

for opt, arg in opts:
    if opt == '-d':
        report_folder = arg


report_files = glob.glob(report_folder + '/ijump_*')    # collect report files
sample_list = list()            # list of samples from reports
summary_table = pd.DataFrame()

# collect all variants
for report in report_files:
    match = re.search('[-|_](\d\w)\.txt', report)
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
    match = re.search('[-|_](\d\w)\.txt', report)
    sample_name = match.group(1)

    for i in range(len(report_df)):
    	if report_df.iloc[i]['Depth'] > 10:
        	ind = str(report_df.iloc[i]['Start']) + report_df.iloc[i]['IS Name'] + str(report_df.iloc[i]['Stop'])
        	summary_table.at[ind, sample_name] = report_df.iloc[i]['Frequency']

summary_table.to_csv('ijump_summary.txt', sep='\t', index=False)

