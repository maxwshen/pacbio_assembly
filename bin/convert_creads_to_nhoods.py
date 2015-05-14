import sys, string, datetime, random, copy, os, commands, fnmatch
import numpy as np
from collections import defaultdict
import read_fasta as rf
import find_read
import itec4

def main():
  print 'Use itec4.py to call this function'
  return


def convert_creads_to_nhoods(reads_file, creads_file, ktmer_headers_file):
  out_file = '/' + '/'.join(creads_file.split('/')[:-1]) + '/nhoods_' + creads_file.split('/')[-1]
  creads = itec4.build_creads_dict(creads_file, reads_file)
  headers = itec4.build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)

  new_text = ''
  for i in range(len(hr)):
    if i % 5000 == 0:
      print i, datetime.datetime.now()
    header = hr[i]
    nh = itec4.get_1_deg_nhood(header, creads, headers)
    new_text += str(i) + ' '
    neighbors_indices = []
    for neighbor_header in nh:
      neighbors_indices.append(hr.index(neighbor_header))
    new_text += ' '.join([str(s) for s in neighbors_indices])
    new_text += '\n'

  with open(out_file, 'w') as f:
    f.write(new_text)

  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
