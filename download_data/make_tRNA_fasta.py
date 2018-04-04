#!/usr/bin/env python

from Bio import SeqIO
import re
from collections import defaultdict
import sys
import six

if len(sys.argv) != 2:
    sys.exit()
ref_path = sys.argv[1]

info_table = ref_path + '/hg38_tRNA.info'
name_dict = {}
with open(info_table, 'r') as info:
    for line in info:
        fields = line.strip().split('\t')
        tRNA_id = re.sub('-[0-9]$','', fields[-2])
        tRNA_name = 'Homo_sapiens_' + fields[-1]
        name_dict[tRNA_name] = tRNA_id
print('made tRNA name hash table', file=sys.stderr)


fasta_file = ref_path + '/hg38-tRNAs.fa'
# remove duplicates and intron
out_seq = defaultdict(set)
out_id = set()
with open(fasta_file, 'r') as fasta:
    for record in SeqIO.parse(fasta, 'fasta'):
        sequence = str(record.seq)
        sequence = re.sub('[actgn]','', sequence)
        name = record.id
        try:
            seq_id = name_dict[name]
        except KeyError:
            seq_id = name_dict['Homo_sapiens_tRNA-Leu-CAA-6-1']
            seq_id = seq_id.replace('-6','-5')
        sequence_set = out_seq[seq_id]
        if sequence not in sequence_set:
            number = len(sequence_set) + 1
            print('>%s-%i\n%sCCAA' %(seq_id, number, sequence))
            out_seq[seq_id].add(sequence)

for key, seqs in six.iteritems(out_seq):
    diff_seq = len(list(seqs))
    if diff_seq > 1:
        print('%s have %i sequences' %(key, diff_seq), file = sys.stderr)
