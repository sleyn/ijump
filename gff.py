import re


class gff:
    def __init__(self, gff_name):
        self.gff_name = gff_name
        self.gff_pos = dict()   # annotation to every position chrom/pos => [lt, name, ID, product]
        self.ann_pos = dict()   # positions for every annotation. Keys are chrom/ids => [ann, contig, start, stop]
        self.capitalize = False

    # read GFF file
    def readgff(self):
        print('Read GFF file ' + self.gff_name)
        gff_file = open(self.gff_name, 'r')
        gff_content = gff_file.read().split('\n##FASTA')    # remove FASTA section
        gff_content = gff_content[0]
        contig_prop_strings = re.findall('##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)', gff_content)  # seqrch for contig properties
        gff_content_sorted = list()

        # sort gff content so headers are grouped with annotations
        for c in contig_prop_strings:
            c_name = re.sub('[^\|]+\|', '', c[0])   # contig name, remove head additional characters separated by "|"
            gff_ges = re.findall('^' + re.escape(c_name) + '\s[^\n]+', gff_content, re.MULTILINE) # take all annotations from the current chromosome
            x = [' '.join(c)]
            x.extend(gff_ges)
            gff_content_sorted.append('\n'.join(x))

        for contig in gff_content_sorted:
            region = re.search('\|?([^\s\|]+)\s+\d+\s+(\d+)', contig)
            chrom = region.group(1)
            if self.capitalize:
                chrom = chrom.upper()

            clen = int(region.group(2))
            self.gff_pos[chrom] = list()
            print('Contig: ' + chrom)

            for pos in range(clen + 1):     # populate all positions with the empty lists
                self.gff_pos[chrom].append(list())

            items = contig.split('\n')

            if len(items) == 1:     # if contig do not have items - fill it with ''
                for i in range(clen + 1):
                    self.gff_pos[chrom][i] = ['', '', '', '']
                continue

            items = [x for x in items if x != '']    # remove blank lines in gff file elements

            # for intergenic space annotations
            prev_pos = 0
            prev_ann = ['', '', '', '']
            prev_orient = '+'

            # for annotation of intergenic space near origin
            first_pos = 0
            last_pos = clen
            first_ann = ''
            last_ann = ''
            first_orient = '+'
            last_orient = '+'

            for j in range(1, len(items)):      # process each annotation
                fields = items[j].split('\t')

                id = '-'   # get id
                id_match = re.search('ID=([^;]+)(;|$)', fields[8])
                if id_match is not None:
                    id = id_match.group(1)

                lt = '-'  # get locus tag
                lt_match = re.search('locus_tag=([^;]+)(;|$)', fields[8])
                if lt_match is not None:
                    lt = lt_match.group(1)

                name = '-'  # get trivial name
                name_match = re.search('gene=([^;]+)(;|$)', fields[8])
                if name_match is not None:
                    name = name_match.group(1)

                product = '-'  # get product
                product_match = re.search('product=([^;]+)(;|$)', fields[8])
                if product_match is not None:
                    product = product_match.group(1)

                # fill positions with information
                for k in range(int(fields[3]), int(fields[4]) + 1):
                    if self.gff_pos[chrom][k]:
                        self.gff_pos[chrom][k][0] += ';' + lt
                        self.gff_pos[chrom][k][1] += ';' + name
                        self.gff_pos[chrom][k][2] += ';' + id
                        self.gff_pos[chrom][k][3] += ';' + product
                    else:
                        self.gff_pos[chrom][k] = [lt, name, id, product]

                # annotate integenic space
                if j > 1 and prev_pos < int(fields[3]):
                    orient1 = 'ds'  # downstream
                    orient2 = 'us'  # upstream

                    if prev_orient is '+':
                        orient1 = 'ds'
                    elif prev_orient is '-':
                        orient1 = 'us'

                    if fields[6] is '-':
                        orient2 = 'ds'
                    elif fields[6] is '+':
                        orient2 = 'us'

                    for k in range(prev_pos + 1, int(fields[3])):
                        if not self.gff_pos[chrom][k]:
                            self.gff_pos[chrom][k] = ['', '', '', '']
                        self.gff_pos[chrom][k][0] = orient1 + '_' + prev_ann[0] + '<>' + orient2 + \
                                                    '_' + self.gff_pos[chrom][int(fields[3])][0]  # for locus_tag
                        self.gff_pos[chrom][k][1] = orient1 + '_' + prev_ann[1] + '<>' + orient2 + \
                                                    '_' + self.gff_pos[chrom][int(fields[3])][1]  # for gene name
                        self.gff_pos[chrom][k][2] = orient1 + '_' + prev_ann[2] + '<>' + orient2 + \
                                                    '_' + self.gff_pos[chrom][int(fields[3])][2]  # for ID
                        self.gff_pos[chrom][k][3] = orient1 + '_' + prev_ann[3] + '<>' + orient2 + \
                                                    '_' + self.gff_pos[chrom][int(fields[3])][3]  # for product

                prev_pos = int(fields[4])
                prev_ann = self.gff_pos[chrom][int(fields[4])]
                prev_orient = fields[6]

                if j == 1:
                    first_pos = int(fields[3])
                    if fields[6] is '+':
                        first_orient = 'us'
                    elif fields[6] is '-':
                        first_orient = 'ds'
                    first_ann = self.gff_pos[chrom][int(fields[3])].copy()

                if j == len(items) - 1:
                    if fields[6] is '-':
                        last_orient = 'us'
                    elif fields[6] is '+':
                        last_orient = 'ds'
                    last_pos = int(fields[4])
                    last_ann = self.gff_pos[chrom][int(fields[4])].copy()

            # annotate between 1 and first gene
            for j in range(first_pos):
                if not self.gff_pos[chrom][j]:
                    self.gff_pos[chrom][j] = ['', '', '', '']
                self.gff_pos[chrom][j][0] = last_orient + '_' + last_ann[0] + '<>' + first_orient + '_' + first_ann[0]
                self.gff_pos[chrom][j][1] = last_orient + '_' + last_ann[1] + '<>' + first_orient + '_' + first_ann[1]
                self.gff_pos[chrom][j][2] = last_orient + '_' + last_ann[2] + '<>' + first_orient + '_' + first_ann[2]
                self.gff_pos[chrom][j][3] = last_orient + '_' + last_ann[3] + '<>' + first_orient + '_' + first_ann[3]

            # annotate between last gene and the end of contig
            for j in range(last_pos, clen + 1):
                if not self.gff_pos[chrom][j]:
                    self.gff_pos[chrom][j] = ['', '', '', '']
                self.gff_pos[chrom][j][0] = last_orient + '_' + last_ann[0] + '<>' + first_orient + '_' + first_ann[0]
                self.gff_pos[chrom][j][1] = last_orient + '_' + last_ann[1] + '<>' + first_orient + '_' + first_ann[1]
                self.gff_pos[chrom][j][2] = last_orient + '_' + last_ann[2] + '<>' + first_orient + '_' + first_ann[2]
                self.gff_pos[chrom][j][3] = last_orient + '_' + last_ann[3] + '<>' + first_orient + '_' + first_ann[3]

    # convert gff_pos => ann to ann => pos for locus tags
    # id => locus tag, contig, start, stop
    def pos_to_ann(self):
        print('Convert GFF table')
        ann_id = 0
        for contig in self.gff_pos.keys():
            prev_ann = ''
            self.gff_pos[contig].append(['-', '-'])                             # add one position to the end for
            self.ann_pos[contig] = dict()
                                                                                # porpose of an algorithm
            for pos in range(len(self.gff_pos[contig]) - 1):
                if prev_ann != self.gff_pos[contig][pos][0]:
                    prev_ann = self.gff_pos[contig][pos][0]
                    if ann_id in self.ann_pos[contig]:
                        self.ann_pos[contig][ann_id][3] = pos
                    ann_id += 1
                    self.ann_pos[contig][ann_id] = [self.gff_pos[contig][pos][0], contig, pos, pos]
