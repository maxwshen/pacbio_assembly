# Combines error corrected read sets into their contigs.
# Input: ec_contigs.fasta
# Format:
#   > (contig number)
#   ACTG... (seq)

# A-Bruijn Graph Construction

import sys, string, datetime, random, copy, os, commands
import numpy as np
from collections import defaultdict
import find_read

def main():
  ec_contig_file = '/home/mshen/research/ec_contigs_900.fasta'
  out_file = '/home/mshen/research/ec_contigs_combined_900best.fasta'
  combine_ec_contigs(ec_contig_file, out_file)

def combine_ec_contigs(ec_contig_file, out_file):
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'
  contigs = []
  with open(ec_contig_file) as f:
    curr_num = -1
    curr_contig = ''
    start_new_contig = False
    build_contig = False    
    for i, line in enumerate(f):
      if start_new_contig:
        if len(curr_contig) > 0:
          contigs.append(curr_contig)
        start_new_contig = False
        build_contig = True
        curr_contig = line.strip()
      if line[0] == '>':
        num = int(line.strip().replace('>', ''))
        if num > curr_num:
          start_new_contig = True
          curr_num = num
          build_contig = False
      if build_contig and line[0] != '>':
        temp_contig_file = 'temp_contig.fasta'
        with open(temp_contig_file, 'w') as f:
          f.write('>' + str(curr_num) + '\n' + curr_contig)
        next_read_file = 'temp_nextread.fasta'
        with open(next_read_file, 'w') as f:
          f.write('>temp\n' + line)
        blasr_out = commands.getstatusoutput(blasr_exe + ' ' + next_read_file + ' ' + temp_contig_file + ' ' + blasr_options)[1]
        print blasr_out
        if len(blasr_out) == 0:
          build_contig = False
        else:
          end = int(blasr_out.split()[10])
          curr_contig += line.strip()[end:]

  with open(out_file, 'w') as f:
    for i in range(len(contigs)):
      f.write('>' + str(i) + '\n' + contigs[i] + '\n')
  return

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start