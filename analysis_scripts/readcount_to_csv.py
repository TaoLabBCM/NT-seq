#!/usr/bin/python

import sys
import re
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd

readcount_fwd = sys.argv[1]
readcount_rev = sys.argv[2]
ref_genome = sys.argv[3]
csv_data_output = sys.argv[4]

def wig_to_datafram(wig_file):
    wig_dict = {}
    wig_dict['chrom'] = []
    wig_dict['pos'] = []
    wig_dict['A'] = []
    wig_dict['C'] = []
    wig_dict['G'] = []
    wig_dict['T'] = []
    with open(wig_file,'r') as wig_f:
        for line in wig_f:
            line.strip()
            if not line:
                continue
            if line.startswith('track'):
                continue
            if line.startswith('#'):
                continue
            if line.startswith('variable'):
                chrom_name = line.split(' ')[1].split('=')[1]
                continue
            wig_dict['chrom'].append(chrom_name)
            wig_dict['pos'].append(line.split('\t')[0])
            wig_dict['A'].append(line.split('\t')[1])
            wig_dict['C'].append(line.split('\t')[2])
            wig_dict['G'].append(line.split('\t')[3])
            wig_dict['T'].append(line.split('\t')[4])
    return(pd.DataFrame(data=wig_dict))
ref_fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(ref_genome, "fasta")}
readcount_fwd_data = wig_to_datafram(readcount_fwd)
readcount_fwd_data['strand'] = '+'
readcount_rev_data = wig_to_datafram(readcount_rev)
readcount_rev_data['strand'] = '-'
csv_data = pd.concat([readcount_fwd_data,readcount_rev_data])
col_list=['A','C','G','T']
csv_data[col_list] = csv_data[col_list].astype(float)
csv_data['cor'] = csv_data[col_list].sum(axis=1)
ref_base = []
for i in range(len(csv_data)):
    chrom_name = csv_data.iloc[i,0]
    base_index = int(csv_data.iloc[i,1])-1
    ref_base.append(str(ref_fasta_dict[chrom_name])[base_index])
csv_data['ref'] = ref_base
csv_data[['chrom','pos','ref','strand','A','C','G','T','cor']].to_csv(csv_data_output,index=False)