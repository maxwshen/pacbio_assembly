# Input: 2 blasr files of reads to the full genome, in m1 format (concise). Should be ~5mb
# Example:
# m120114_011938_42177_c100247042550000001523002504251220_s1_p0/8/0_1360/0_1360/0_1296 gi|556503834|ref|NC_000913.3| 0 0 -6256 99.9202 2822360 2823612 4641652 42 1295 1296 1024

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf
import find_read
import kmer_matching
import generate_new_reads
from collections import defaultdict

def main():
  base = '/home/mshen/research/data/READS_blasr_iter_0.out'
  ecf = sys.argv[1]

  quantify_ec(base, ecf)

def quantify_ec(base, ecf):
  base_data = defaultdict(list)  # Key = header, Value = [accuracy, start_pos, read_length]
  with open(base) as f:
    for i, line in enumerate(f):
      words = line.split()
      header = '/'.join(words[0].split('/')[:-1])
      accuracy = float(words[5])
      start_pos = int(words[6])
      read_len = int(words[11])
      base_data[header] = [accuracy, start_pos, read_len]

  unfound_headers = []
  new_data = defaultdict(list)
  with open(ecf) as f:
    for i, line in enumerate(f):
      if i % 10000 == 0:
        print i, datetime.datetime.now()
      words = line.split()
      header = '/'.join(words[0].split('/')[:-1])
      if header in base_data.keys():
        accuracy = float(words[5]) - base_data[header][0]
        start_pos_diff = int(words[6]) - base_data[header][1]
        read_len_diff = int(words[11]) - base_data[header][2]
        read_len = int(words[11])
        new_data[header] = [accuracy, start_pos_diff, read_len_diff, read_len]
      else:
        unfound_headers.append(header)

  better_acc = []
  num_unchanged_acc = 0
  worse_acc = []
  start_pos_bigdiff = []
  read_len_pctdiff = []
  for k in new_data.keys():
    # Accuracy
    if new_data[k][0] > 0:
      better_acc.append(new_data[k][0])
    if new_data[k][0] == 0:
      num_unchanged_acc += 1
    if new_data[k][0] < 0:
      worse_acc.append(new_data[k][0])

    # Start position difference
    if new_data[k][1] > 10000 or new_data[k][1] < -10000:
      start_pos_bigdiff.append(new_data[k][1])

    # Read length percent difference
    if new_data[k][2] < 0:
      read_len_pctdiff.append(float(new_data[k][2]) / float(new_data[k][3] - new_data[k][2]))
    else:
      read_len_pctdiff.append(float(new_data[k][2]) / float(new_data[k][3] + new_data[k][2]))

  print 'Unfound headers:', len(unfound_headers)
  print 'Better accuracy:', np.average(better_acc), np.std(better_acc), len(better_acc)
  print 'Unchanged accuracy:', num_unchanged_acc
  print 'Decreased accuracy:', np.average(worse_acc), np.std(worse_acc), len(worse_acc)
  print 'Changed position:', len(start_pos_bigdiff) 
  print 'Avg. Readlen Pct. Change:', np.average(read_len_pctdiff) * 100
  return


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start