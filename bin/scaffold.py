# Scaffolds together combined contigs in the same folder

import sys, string, datetime, random, copy, os, commands, fnmatch, re
import numpy as np
from collections import defaultdict
import itec4
import mylib

def main():
  contigs_dr = sys.argv[1]
  
  scaffold_main(contigs_dr)

  return
  

def scaffold_main(contigs_dr):
  cns = []
  num_contigs = 0
  for fn in os.listdir(contigs_dr):
    if fnmatch.fnmatch(fn, '*_combined.fasta'):
      cns.append(Contig(fn, num_contigs))
      num_contigs += 1

  for jump in range(1, len(cns) - 1):
    for i in range(0, len(cns) - jump):
      print i, i + jump
      res = attempt_scaffold(cns[i].ref, cns[i + jump].ref, contigs_dr)
      if len(res) != 0:
        print 'success', res
        cns[i].ref = res
        cns[i + jump].ref = res
      new_tried = cns[i].tried + cns[i + jump].tried
      for nt in new_tried:
        cns[nt].tried = new_tried

def attempt_scaffold(c1_fn, c2_fn, contigs_dr):
  # Returns '' upon failure
  # Returns new scaffold name upon success
  acc_cutoff = 0
  dist_to_end = 2000

  blasr_exe = 'blasr'
  blasr_zero = 4      # 0 on debruijn, 4 on Yu's computer
  blasr_zero_len = 8  # 0 on debruijn, 8 on Yu's computer
  blasr_options = '-bestn 1 -m 1 -maxMatch 20'

  c1h, c1r = mylib.read_fasta(contigs_dr + c1_fn)
  c2h, c2r = mylib.read_fasta(contigs_dr + c2_fn)
  c1r = c1r[0]
  c2r = c2r[0]
  
  print 'test'
  status = commands.getstatusoutput(blasr_exe + ' ' + contigs_dr + c1_fn + ' ' + contigs_dr + c2_fn + ' ' + blasr_options)[1]
  print status
  if len(status.split()) == blasr_zero_len:
    print 'FAILED BLASR ALIGNMENT'
    return ''
  else:
    print status                      # TESTING
    acc = float(status.split()[blasr_zero + 5])
    beg_align_r1 = int(status.split()[blasr_zero + 6])
    end_align_r1 = int(status.split()[blasr_zero + 7])
    total_len_r1 = int(status.split()[blasr_zero + 8])
    end_pos_r1 = total_len_r1 - end_align_r1
    beg_align_r2 = int(status.split()[blasr_zero + 9])
    end_align_r2 = int(status.split()[blasr_zero + 10])
    total_len_r2 = int(status.split()[blasr_zero + 11])
    end_pos_r2 = total_len_r2 - end_align_r2

    new_base = ''
    if end_pos_r1 < dist_to_end: # and beg_align_r2 < dist_to_end:
      if acc > acc_cutoff:
        new_base = c1r[: end_align_r1]
        new_base += c2r[end_align_r2 : ]
    elif end_pos_r2 < dist_to_end:
      if acc > acc_cutoff:
        new_base = c2r[: end_align_r2]
        new_base += c1r[end_align_r1 : ]
    if len(new_base) != 0:
      c1n = c1_fn.split('_')[1]
      c2n = c2_fn.split('_')[1]
      if c1_fn.split('_')[-1] == '_scaffold.fasta':
        commands.getstatusoutput('rm -rf ' + c1_fn)
      if c2_fn.split('_')[-1] == '_scaffold.fasta':
        commands.getstatusoutput('rm -rf ' + c2_fn)
      newc_fn = 'contig_' + c1n + c2n + '_scaffold.fasta'
      with open(contigs_dr + newc_fn, 'w') as f:
        f.write('>scaffold\n' + newc_fn)
      return newc_fn
    return ''

def check_done(cns):
  for cn in cns:
    if not cn.is_complete(len(cns)):
      return False
  return True

class Contig():
  def __init__(self, fn, num):
    self.num = num
    self.tried = [num]
    self.neighbors = [num]
    self.ref = fn

  def is_complete(self, cap):
    return len(self.tried) == cap


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
