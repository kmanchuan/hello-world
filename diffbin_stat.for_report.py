#!/usr/bin/env python3

import os
import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Spot:
    def __init__(self, matrix_line=None, x=None, y=None):
        if matrix_line:
            lst = matrix_line.strip('\n').split('\t')
            gene_id, x, y, umi_count, exon_count = lst
            self.genes = [gene_id]
            self.umi_count = int(umi_count)
        elif not matrix_line:
            self.genes = []
            self.umi_count = 0
        self.x = int(x)
        self.y = int(y)

class Matrix(list):
    def __init__(self, arg=None):
        if isinstance(arg, list):
            self.extend(arg)
        elif isinstance(arg, str) and os.path.exists(arg):
            if arg.endswith('.gz'):
                r = gzip.open(arg, 'rt')
            else:
                r = open(arg)
            for line in r:
                if line.startswith(('#', 'geneID')):
                    continue
                rec = Spot(matrix_line=line)
                self.append(rec)
            r.close()
    def binning(self, bin_size=1):
        tmpdict = dict()
        for rec in self:
            x = rec.x // bin_size
            if rec.x % bin_size > 0:
                x += 1
            y = rec.y // bin_size
            if rec.y % bin_size > 0:
                y += 1
            coord = (x, y)
            if not tmpdict.get(coord):
                tmpdict[coord] = Spot(x=x, y=y)
            tmpdict[coord].genes.extend(rec.genes)
            tmpdict[coord].umi_count += rec.umi_count
        grouped_matrix = Matrix(list(tmpdict.values()))
        return grouped_matrix
    def stat(self):
        gene_counts = [len(set(x.genes)) for x in self]
        mean_gene = np.mean(gene_counts)
        median_gene = np.median(gene_counts)
        umi_counts = [x.umi_count for x in self]
        mean_umi = np.mean(umi_counts)
        median_umi = np.median(umi_counts)
        return mean_gene, median_gene, mean_umi, median_umi, gene_counts, umi_counts

if __name__ == '__main__':

    gene_matrix = Matrix(arg=sys.argv[1])
    SN = sys.argv[2]
    #prefix = os.path.splitext(os.path.basename(sys.argv[1]))[0]

    binSizeList = [20, 50, 100, 150, 200]
#    header = ['binSize', 'Mean Gene Types', 'Median Gene Types', 'Mean MIDs', 'Median MIDs']
#    sys.stdout.write('\t'.join(header) + '\n')
    sys.stdout.write(SN)
    fig_size = []
    fig_gene = []
    fig_umi = []
    for size in binSizeList:
        if size == 1:
            stat = gene_matrix.stat()
        else:
            grouped_matrix = gene_matrix.binning(bin_size=size)
            stat = grouped_matrix.stat()
        gene_counts, umi_counts = stat[-2:]
        line = '\t{}\t'.format(size)
        #line += '\t'.join(['{:.2f}'.format(x) for x in stat[:-2]])
        line += '\t'.join(['{:,}'.format(int(x)) for x in stat[:-2]])
        sys.stdout.write(line + '\n')
        fig_size.extend([size] * len(gene_counts))
        fig_gene.extend(gene_counts)
        fig_umi.extend(umi_counts)

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(9, 4), sharex=True)
    sns.violinplot(x=fig_size, y=fig_gene, ax=ax1, scale='width')
    ax1.set_ylabel('gene')
    sns.violinplot(x=fig_size, y=fig_umi, ax=ax2, scale='width')
    ax2.set_ylabel('MID')
    ax2.set_xlabel('binSize')
    fig.savefig('{}_vln.png'.format(SN), dpi=300)


