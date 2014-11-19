import sys
import string
import datetime
import random
import copy
import os
import commands

from collections import defaultdict

def main():
  # input_folder = '/home/lin/nhood_files/'
  input_folder = '/home/mshen/research/yu_ec_22.4_500_nhoods/'

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  # blasr_options = '-bestn 1 -m 1'   # Concise output
  blasr_options = '-bestn 1 -m 0'   # Verbose output
  for name in os.listdir(input_folder):
    print commands.getstatusoutput(blasr_exe + ' ' + input_folder + '/' + name + ' ' + e_coli_genome + ' ' + blasr_options)[1]
    # out = commands.getstatusoutput(blasr_exe + ' ' + input_folder + '/' + name + ' ' + e_coli_genome + ' ' + blasr_options)[1]
    # print out.split()[4] + '\n'
  return

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start