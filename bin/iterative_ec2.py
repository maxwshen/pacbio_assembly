# Performs iterative error correction of reads.
# Beginning with an arbitrary neighborhood,
#   1. Get the full reads from the neighborhood
#   2. Error correct all of these reads, producing a consensus
#   3. Find the farthest kt-mer (with an nhood) in the consensus and repeat
#     - If the coverage is too low or too high, then use k-mer matching instead
# Choose a new nhood when k-mer matching fails to find sufficient reads

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf
import find_read
import kmer_matching
import generate_new_reads
from collections import defaultdict

def main():
  # reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  reads_file = '/home/mshen/research/extracts/extracted_reads_c1050000_s100000.fasta'
  nhood_fold = '/home/mshen/research/e_coli_nh_22.4_unwhole_100k/'

  iterative_ec(reads_file, nhood_fold)

def iterative_ec(reads_file, nhood_fold):
  reads_h, reads_r = rf.read_fasta(reads_file)
  ec_tool = '/home/mshen/research/bin/error_correction.sh'
  # e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  e_coli_genome = '/home/mshen/research/extracts/extracted_genome_c1050000_s100000.fasta'
  _k = 15
  cutoff = 5
  
  ktmers = build_ktmer_dict(nhood_fold)

  traversed = set()
  nhoods = os.listdir(nhood_fold)
  random.shuffle(nhoods)
  high_cutoff = 50
  low_cutoff = 2

  curr_name = ''
  while len(traversed) < len(ktmers):
    if curr_name == '':
      curr_name = random.choice(nhoods)
      while curr_name in traversed:
        curr_name = random.choice(nhoods)

    traversed.add(curr_name)
    curr = nhood_fold + '/' + curr_name
    with open(curr) as f:
      num_lines = len(f.readlines())
    temp_nhreads_f = 'temp.fasta'
    print len(traversed), num_lines
    if low_cutoff < num_lines < high_cutoff:
      get_reads_from_nhood(curr, reads_file, temp_nhreads_f)
    else:
      fr_file = 'temp_fr.fasta'
      get_first_read(curr, reads_file, fr_file)
      _k = 15
      km_cutoff = 10
      kmer_matching.kmer_matching(fr_file, reads_file, _k, km_cutoff, temp_nhreads_f)
    ec_out = commands.getstatusoutput(ec_tool + ' ' + temp_nhreads_f)[1]
    if 'ERROR' in ec_out:
      print 'COULD NOT ERROR CORRECT', ec_out
      curr_name = ''
      continue

    rr_outf = 'temp_replacereads.txt'
    consensus_file = 'corrected_' + temp_nhreads_f
    replace_reads_blasr(consensus_file, temp_nhreads_f, rr_outf)
    generate_new_reads(rr_outf, reads_h, reads_r)

    ec_h, ec_r = rf.read_fasta(consensus_file)

    best_name = ''
    largest_index = 0
    for name in nhoods:
      if name not in traversed:
        ktmer = name.split('_')[2].split('.')[0]
        if ktmer in ec_r[0]:
          if ec_r[0].index(ktmer) > largest_index:
            largest_index = ec_r[0].index(ktmer)
            best_name = name
            print best_name, largest_index
    curr_name = best_name

def get_first_read(fn, reads_file, out_file):
  with open(fn) as f:
    h = f.readlines()[0].strip()
  with open(out_file, 'w+') as f2:
    f2.write(find_read.find_read(h, reads_file))
  return

def get_reads_from_nhood(fn, reads_file, out_file):
  # Takes about 1 second for nhoods with 20 reads
  h = []
  with open(fn) as f:
    with open(out_file, 'w+') as f2:
      for i, line in enumerate(f):
        if line[0] == '>':
          f2.write(find_read.find_read(line.strip(), reads_file))
  return

def build_ktmer_dict(nhood_fold):
  # Template: nhood_nh1091437_ATGATCAACCTGAATACGCTGG.fasta
  ktmers = dict()   # Key = ktmer, dict = file name
  for name in os.listdir(nhood_fold):
    ktmer = name.split('_')[2].split('.')[0]
    ktmers[ktmer] = name
  return ktmers


def replace_reads_blasr(consensus_file, km_file_name, out_file):
  # Given a set of error corrected consensuses derived from corresponding nhoods
  # (they share a name), use blasr to align the EC consensus to its nhood and output
  # in text format the indices that should be used for replacement.

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'
  accuracy_cutoff = 80
  corrections = []
  heads, reads = rf.read_fasta(km_file_name)

  for i in range(len(reads)):
    h = heads[i]
    r = reads[i]
    with open('temp.txt', 'w+') as f:
      f.write(h + '\n' + r)
    blasr_out = commands.getstatusoutput(blasr_exe + ' ' + consensus_file + ' ' + 'temp.txt' + ' ' + blasr_options)[1]
    if len(blasr_out) > 0:
      accuracy = blasr_out.split()[5]
      nh_beg = blasr_out.split()[6]
      nh_end = blasr_out.split()[7]
      ec_beg = blasr_out.split()[9]
      ec_end = blasr_out.split()[10]
      if float(accuracy) > accuracy_cutoff:
        info = h + ' ' + nh_beg + ' ' + nh_end + ' ' + consensus_file + ' ' + ec_beg + ' ' + ec_end
        corrections.append(info)
  with open(out_file, 'w+') as f:
    f.write('\n'.join(corrections))
  return


def generate_new_reads(rr_file, reads_h, reads_r):
  # Replaces EC consensus into reads
  lib = defaultdict(list)
  with open(rr_file) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        h = line.split()[0]
        lib[h].append(' '.join(line.split()[1:]))

  for i in range(len(reads_r)):
    if reads_h[i] in lib.keys():
      readl = list(reads_r[i])
      for info in lib[h]:
        sinfo = info.split()
        read_beg = int(sinfo[0])
        read_end = int(sinfo[1])
        ec_file = sinfo[2]
        ec_beg = int(sinfo[3])
        ec_end = int(sinfo[4])
        with open(ec_file) as ec_f:
          ec_read = ec_f.readlines()[1].strip()
        readl[read_beg : read_end + 1] = ec_read[ec_beg : ec_end + 1]
      reads_r[i] = ''.join(readl)
  return


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start