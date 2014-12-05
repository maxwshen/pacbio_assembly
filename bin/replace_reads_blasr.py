# replace_reads.py is a program that finds where in each read (header) a nhood should be replaced in
#
# Output Format:
# >m120114_011938_42177_c100247042550000001523002504251220_s1_p0/71370/0_3131/0_3131 1643 2157 /home/mshen/research/yu_ec_22.4_500_nhoods/nhood_nh2943330_AGTACCATGAACGTTTTTAATCC.fasta_Cov5.fasta 9 485
# Header, read_start_pos, read_end_pos, ec_nhood_file_name, ec_start, ec_end 

import sys, string, datetime,random, copy, os, commands
import find_read

from collections import defaultdict

def main():
  # input_folder = '/home/lin/nhood_files/'
  nhoods_fold = '/home/mshen/research/e_coli_nhoods_500_22.4_unwhole/'
  ec_fold = '/home/mshen/research/yu_ec_22.4_500_nhoods/'

  replace_reads_blasr(nhoods_fold, ec_fold)
  return


def replace_reads_blasr(nhoods_fold, ec_fold, folder = False):
  # nhoods_fold not needed if we're using k-mer matching
  # Given a set of error corrected consensuses derived from corresponding nhoods
  # (they share a name), use blasr to align the EC consensus to its nhood and output
  # in text format the indices that should be used for replacement.

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  km_file_name = '/home/mshen/research/nohup_km_22.4_fullgenome.out'

  ec_names = os.listdir(ec_fold)
  nh_names = os.listdir(nhoods_fold)
  d = build_hash(nh_names, ec_names)

  accuracy_cutoff = 80

  corrections = []

  for key in d:
    # headers = get_headers_nhood(nhoods_fold + d[key])
    headers = get_headers_km(km_file_name, key)
    for h in headers:
      full_read = find_read.find_read(h, reads_file)
      with open('temp.txt', 'w+') as f:
        f.write(full_read)
      blasr_out = commands.getstatusoutput(blasr_exe + ' ' + ec_fold + key + ' ' + 'temp.txt' + ' ' + blasr_options)[1]
      if len(blasr_out) > 0:
        accuracy = blasr_out.split()[5]
        nh_beg = blasr_out.split()[6]
        nh_end = blasr_out.split()[7]
        ec_beg = blasr_out.split()[9]
        ec_end = blasr_out.split()[10]
        if float(accuracy) > accuracy_cutoff:
          info = h + ' ' + nh_beg + ' ' + nh_end + ' ' + ec_fold + key + ' ' + ec_beg + ' ' + ec_end
          print info
          corrections.append(info)
  return corrections

def get_headers_km(km_file_name, ec_name):
  # Inefficient, opens file every time. But structurally similar to get_headers_nhood
  get_headers = False
  headers = []
  with open(km_file_name) as f:
    for i, line in enumerate(f):
      if get_headers:
        if line[0] == '>':
          headers.append(line.strip())
        if len(line.strip()) == 0:
          break
      if line.split('/')[-1].strip() == ec_name:
        get_headers = True
  return headers

def get_headers_nhood(file_name):
  headers = []
  with open(file_name) as f:
    lines = f.readlines()
    for line in lines:
      if line[0] == '>':
        headers.append(line.strip())
  return headers

def build_hash(nh_names, ec_names):
  # Builds a dict. Key = ec_file_name. Value = nh_name
  nh_names = sorted(nh_names)
  ec_names = sorted(ec_names)
  found = False
  d = dict()
  for i in range(len(ec_names)):
    for nh_name in nh_names:
      if ec_names[i][:20] == nh_name[:20]:
        d[ec_names[i]] = nh_name
        found = True
        break
    if found:
      nh_names.remove(d[ec_names[i]])
      found = False
  return d


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start