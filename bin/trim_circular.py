# Trims a combined contig to make it circular

import sys, string, datetime, random, copy, os, commands, fnmatch, re
import numpy as np
from collections import defaultdict
import mylib

blasr_exe = 'blasr'
blasr_zero = 4      # 0 on debruijn, 4 on Yu's computer
blasr_zero_len = 8  # 0 on debruijn, 8 on Yu's computer
blasr_options = '-bestn 1 -m 1 -maxMatch 20'   # Concise output

def main():
  cc_fn = sys.argv[1]
  res = trim_circular(cc_fn)
  while res != '':
    next_res = trim_circular(cc_fn)
    status = commands.getstatusoutput('rm -rf ' + res)[1]
    status = commands.getstatusoutput('mv ' + next_res + ' ' res)[1]
    res = next_res
  return

def trim_circular(cc_fn):
  h, r = mylib.read_fasta(cc_fn)
  r = r[0]
  LENALIGN = 50000
  DIST_TO_END = 100
  ACC_CUTOFF = 99.3

  if len(r) < LENALIGN * 2:
    print 'contig shorter than', LENALIGN * 2
    return

  temp1_fn = cc_fn + 'temp1.fasta'
  temp2_fn = cc_fn + 'temp2.fasta'

  with open(temp1_fn, 'w') as f:
    f.write('>temp1\n' + r[-LENALIGN : ])

  with open(temp2_fn, 'w') as f:
    f.write('>temp2\n' + r[:LENALIGN])

  status = commands.getstatusoutput(blasr_exe + ' ' + temp1_fn +' ' + temp2_fn + ' ' + blasr_options)[1]
  if len(status.split()) == blasr_zero_len:
    return False
  print status                    # TESTING
  r2_strand_dir = int(status.split()[blasr_zero + 2])
  r1_strand_dir = int(status.split()[blasr_zero + 3])
  accuracy = float(status.split()[blasr_zero + 5])
  beg_align_r1 = int(status.split()[blasr_zero + 6])
  end_align_r1 = int(status.split()[blasr_zero + 7])
  total_len_r1 = int(status.split()[blasr_zero + 8])
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[blasr_zero + 9])
  end_align_r2 = int(status.split()[blasr_zero + 10])
  total_len_r2 = int(status.split()[blasr_zero + 11])
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2   # Average alignment length

  status = commands.getstatusoutput('rm -rf ' + temp1_fn)[1]
  status = commands.getstatusoutput('rm -rf ' + temp2_fn)[1]

  if beg_align_r1 < DIST_TO_END and end_pos_r2 < DIST_TO_END and accuracy > ACC_CUTOFF:
    new_r = r[:-beg_align_r2]
    out_fn = cc_fn.split('.fasta')[0] + '.circle.fasta'
    with open(out_fn, 'w') as f:
      f.write(h[0] + '\n' + new_r)
    print 'SUCCESS'
    return out_fn
  else:
    print 'FAILURE'
    return ''


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start