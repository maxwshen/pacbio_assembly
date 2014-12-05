# Performs iterative error correction of reads.
# Beginning with an arbitrary error corrected consensus,
#   1. Use k-mer matching to find corresponding reads
#   2. Use blasr to align the EC consensus to these reads, finding indices for replacement
#   3. Replace EC consensus into reads, then grab next 500bp and perform error correction
# This process is repeated until the coverage falls too low for proper error correction

import sys, string, datetime, random, copy, os, commands
import numpy as np

import read_fasta as rf
import find_read
import kmer_matching
import generate_new_reads
from collections import defaultdict

def main():
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  consensus_file = sys.argv[1]


  iterative_ec(consensus_file, reads_file)

def iterative_ec(consensus_file, reads_file):
  reads_h, reads_r = rf.read_fasta(reads_file)
  ec_tool = '/home/mshen/research/bin/error_correction.sh'
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1'
  _k = 15
  cutoff = 10
  extend_len = 1000
  
  build_km_lib(reads_h, reads_r, _k)
  return

  while True:
    km_outf = 'km_' + consensus_file.translate(None, '/') + '.fasta'
    kmer_matching.kmer_matching(consensus_file, reads_file, 15, 5, km_outf)
    print 'kmer matching done'

    rr_outf = 'temp_replacereads.txt'
    replace_reads_blasr(consensus_file, km_outf, rr_outf)
    print 'replace reads calculated'

    next_area_f = 'temp_nextarea.fasta'
    generate_new_reads(rr_outf, reads_h, reads_r, extend_len, next_area_f)
    print 'found next area'

    status = commands.getstatusoutput(ec_tool + ' ' + next_area_f)[1]
    print status
    if 'ERROR' in status:
      sys.exit(0)
    consensus_file = 'corrected_temp_nextarea.fasta'    # Based on Yu's ec_tool
    blasr_out = commands.getstatusoutput(blasr_exe + ' ' + consensus_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]
    print blasr_out

def build_km_lib(reads_h, reads_r, _k):
  # tooo slowwww. probably due to the large amount of memory necessary
  lib = defaultdict(list)   # Key = k-mer, value = list of [header, number of appearances] 
  for i in range(len(reads_r)):
    print i
    read = reads_r[i]
    head = reads_h[i]
    for j in range(len(read) - _k + 1):
      kmer = read[j : j + _k]
      if kmer not in lib.keys():
        lib[kmer].append([head, 1])
      elif head in [s[0] for s in lib[kmer]]:
        ind = [s[0] for s in lib[kmer]].index(head)
        lib[kmer][ind][1] += 1
  print len(lib)
  return lib



def kmer_matching(ec_seq_file, read_file, _k, cutoff, out_file):
  # ec_seq_file should be a fasta file with only one sequence

  hs, rs = rf.read_fasta(ec_seq_file)
  ec_seq = rs[0]

  kmers = set()
  for i in range(len(ec_seq) - _k + 1):
    kmers.add(ec_seq[i:i + _k])

  reads = dict()    # Key = header, value = num shared kmers
  hr, rr = rf.read_fasta(read_file)
  for i in range(len(rr)):
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  headers = []
  for key in sorted(reads, key = reads.get, reverse = True):
    if reads[key] < cutoff:
      break
    headers.append(key)

  print '\n' + ec_seq_file + '\n' + str(len(headers))
  get_reads_from_headers(headers, read_file, out_file)
  # return
  

  # Use blasr to check accuracy of kmer-matching found reads

  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  blasr_options = '-bestn 1'
  blasr_out = commands.getstatusoutput(blasr_exe + ' ' + out_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]

  to_print = []
  for line in blasr_out.splitlines():
    head = '>' + '/'.join(line.split()[0].split('/')[:-1])
    start_pos = line.split()[6]
    end_pos = line.split()[7]
    length = int(end_pos) - int(start_pos)
    info = (head, str(reads[head]), start_pos, end_pos, str(length))
    to_print.append(info)

  for t in sorted(to_print, key = lambda tup: int(tup[1]), reverse = True):
    print t[0]
    print '\t'.join(t[1:])

  return


def generate_new_reads(rr_file, reads_h, reads_r, extend_len, next_area_f):
  # Replaces EC consensus into reads and also produces next 500 or whatever extend_len is
  next_area = []

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
        na = readl[read_end : read_end + extend_len + 1]
        if len(na) > extend_len:
          next_area.append(reads_h[i] + '\n' + ''.join(na))
      reads_r[i] = ''.join(readl)

  with open(next_area_f, 'w+') as f:
    f.write('\n'.join(next_area))

  return


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


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start