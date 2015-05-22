# This finds the percent of (k,t)-mers that are in the genome.
# Checks the reverse complement of the genome

import read_fasta as rf
from collections import defaultdict

def rc(kmer):
  match = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
  return ''.join([match[s] for s in kmer[::-1]])

ktmers_file = '/home/mshen/research/data/20k_v2/temp_ktmer_headersrx_27_6_rc_v2.out'
genome_file = '/home/mshen/research/data/ecoli_consensus_mark.fasta'
# genome_file = '/home/mshen/research/data/e_coli_genome.fasta'


gh, gr = rf.read_fasta(genome_file)
gr = gr[0]
print len(gr)

ktmers_dict = defaultdict(list)
with open(ktmers_file) as f:
  for i, line in enumerate(f):
    ktmer = line.split()[0]
    headers = line.split()[1:]
    ktmers_dict[ktmer] = headers

_k = 27
total = len(ktmers_dict.keys())
num_found = 0
found = set()
for i in range(len(gr) - _k + 1):
  curr_ktmer = gr[i : i + _k]
  if curr_ktmer in ktmers_dict and curr_ktmer not in found:
    num_found += 1
    found.add(curr_ktmer)
  rck = rc(curr_ktmer)
  if rck in ktmers_dict and rck not in found:
    num_found += 1
    found.add(rck)

print num_found, total, float(num_found) / float(total)

