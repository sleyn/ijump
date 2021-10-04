import sys
import getopt
import re
import pandas as pd

html_name = ""

# read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:")
except getopt.GetoptError:
    print("isfinder_parse.py -i <ISfinder BLAST HTML page>")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-i':
        isfinder_name = arg

# check if all argumants were provided

if not isfinder_name:
    print("ISfinder BLAST HTML file is not specified")
    sys.exit(1)


# read ISfinder output html file
isfinder_file = open(isfinder_name, 'r')
isfinder_content = isfinder_file.read()
isfinder_file.close()

# collect information from the table

isfinder_content = re.sub('</article>.+', '', isfinder_content, re.DOTALL) #remove HTML code after BLAST output

isels_final = pd.DataFrame(columns=['contig', 'name', 'family', 'group', 'origin', 'score', 'e-value', 'start', 'stop',
                              'check'])  # specific inforamtion about every hit

# Check if the Query= string exists
query_exists = False

for contig in isfinder_content.split('<b>Query=</b>')[1:]:
    query_exists = True
    isels_ge = pd.DataFrame(columns=['name', 'family', 'group', 'origin'])  # general information about IS
    isels = pd.DataFrame(columns=['contig', 'name', 'family', 'group', 'origin', 'score', 'e-value', 'start', 'stop',
                                  'check'])  # specific inforamtion about every hit

    gname_match = re.match(' (\S+)[^\n]*\n\nLength=(\d+)\n',contig)
    gname = gname_match.group(1)        # contig name
    gsize = int(gname_match.group(2))   #contig size

    for isel in re.findall('<a[^>]+>([^<]+)</a></td><td>([^<]*)</td><td>([^<]*)</td><td><a[^>]+>([^<]*)</a></td><td><a[^>]+>[^<]+</a></td><td>[^<]+</td></tr>',contig):     # IS name, family, group, genome
        isel_df = pd.DataFrame([list(isel)], columns=['name', 'family', 'group', 'origin'])
        isels_ge = isels_ge.append(isel_df, sort=False)

    isels_ge = isels_ge.set_index('name')

    # collect coordinates of IS elements
    match = re.search('<pre>&gt;(.+)', contig, re.DOTALL)
    if match is None:
        continue

    alignments_block = match.group(1)
    for alignments in alignments_block.split('&gt;'):
        match = re.search('</a> (\S+) <span', alignments)
        name = match.group(1)
        index = 1

        for align in alignments.split(' Score')[1:]:
            evalue = float(re.search('Expect = ([^\n]+)', align).group(1))
            if evalue > 1E-30:
                continue
            bit = int(re.search(' = (\d+) bits', align).group(1))
            match = re.findall('Query[^\n]+\n',align)
            smatch = re.search('Query\s+(\d+)\s+',match[0])
            start = smatch.group(1)
            smatch = re.search('Query\s+\d+\s+[ATGCNatgcn-]+\s+(\d+)', match[-1])
            try:
                stop = smatch.group(1)
            except:
                print(name)
                print(match)

            add_list = [gname, name + '_' + str(index)]
            add_list.extend(list(isels_ge.loc[name,:]))
            add_list.extend([bit, evalue, start, stop, 0])
            isel_df = pd.DataFrame([add_list], columns=['contig', 'name', 'family', 'group', 'origin', 'score', 'e-value', 'start', 'stop', 'check'])
            isels = isels.append(isel_df, sort=False)
            index = index + 1

    isels = isels.sort_values(by=['score'], ascending=False)
    isels = isels.set_index('name')

    genome = [0]*gsize       # set genome coverage list

    for index in isels.index:
        genome_temp = genome.copy()
        is_new = 0                  # new space occupied by IS element
        stop = int(isels.at[index, 'stop'])
        start = int(isels.at[index,'start'])
        is_length = stop - start             # length of the IS element
        for pos in range(start - 1, stop):
            if genome_temp[pos] == 0:
                is_new = is_new + 1

            genome_temp[pos] = 1

        if is_new / is_length > 0.75:           # check if IS element occupies a new region at least 75% of its length
            genome = genome_temp
            isels.at[index, 'check'] = 1

    isels = isels[isels['check']==1]            # remove overlapping IS hits with less score
    isels = isels.drop(columns=['check'])       # remove check column
    isels_final = isels_final.append(isels, sort=False)

if query_exists:
    # write full table in file
    isels_final.to_csv('ISTable_full.txt', sep='\t')
    # write table for further processing to file
    isels_final.to_csv('ISTable_processing.txt', sep='\t', columns=['contig', 'start', 'stop'], header=False)
else:
    print('\nWARNING!:'
          'No "Query=" strings with contig names found in the provided page.\n'
          'ISFinder BLAST requires short Fasta headers to include contig names.\n'
          'Please make shorter contig names or remove auxiliary information from the Fasta header.')
