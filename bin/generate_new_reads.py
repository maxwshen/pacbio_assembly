# Input Format (from replace_reads_blasr.py):
# >m120114_011938_42177_c100247042550000001523002504251220_s1_p0/71370/0_3131/0_3131 1643 2157 /home/mshen/research/yu_ec_22.4_500_nhoods/nhood_nh2943330_AGTACCATGAACGTTTTTAATCC.fasta_Cov5.fasta 9 485
# Header, read_start_pos, read_end_pos, ec_nhood_file_name, ec_start, ec_end 

import sys, string, datetime, random, copy, os, commands
import find_read

from collections import defaultdict

def main():
  replace_reads_output = '/home/mshen/research/nohup_replacereads_22.4_fullgenome.out'
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  new_reads_file = 'NEXT500_3276526_22.4.fasta'

  rr_km_file = '/home/mshen/research/nohup_rr_km_head_22.4.out'
  get_next_500(rr_km_file, reads_file, new_reads_file)
  # generate_new_reads(replace_reads_output, reads_file, new_reads_file)
  return

def get_next_500(rr_km_file, reads_file, new_reads_file):
  # An error-corrected 500bp nhood is aligned to kmer-matching reads
  # We then grab the next 500bp from these reads and output them
  threshold = 450
  dist = 500
  reads = dict()   # Key = header, value = (beg, end)

  with open(rr_km_file) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        words = line.split()
        h = words[0]
        read_beg = int(words[1])
        read_end = int(words[2])
        ec_file = words[3]
        ec_beg = int(words[4])
        ec_end = int(words[5])
        if ec_end - ec_beg > threshold:
          reads[h] = (read_beg, read_end)

  grab_read = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if grab_read:
        grab_read = False
        print h + '\n' + line[end : end + 500].strip()
      if line[0] == '>' and line.strip() in reads.keys():
        grab_read = True
        h = line.strip()
        beg = reads[h][0]
        end = reads[h][1]


def generate_new_reads(rr_file, reads_file, new_reads_file):
  lib = defaultdict(list)
  with open(rr_file) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        h = line.split()[0]
        lib[h].append(' '.join(line.split()[1:]))

  with open(new_reads_file, 'w+') as f:
    pass

  found_header = False
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if found_header:
        found_header = False
        read = line.strip()
        readl = list(read)
        log = dict()    # Key = pos, value = offset
        for info in lib[h]:
          sinfo = info.split()
          read_beg = int(sinfo[0])
          read_end = int(sinfo[1])
          ec_file = sinfo[2]
          ec_beg = int(sinfo[3])
          ec_end = int(sinfo[4])
          with open(ec_file) as ec_f:
            ec_read = ec_f.readlines()[1].strip()
          offset = 0
          for key in log.keys():
            if read_beg > key:
              offset += log[key]
          readl[read_beg - offset: read_end + 1 - offset] = ec_read[ec_beg : ec_end + 1]
          log[read_beg + ec_end - ec_beg] = read_end - read_beg - ec_end + ec_beg
        with open(new_reads_file, 'a') as f:
          f.write(h + '\n' + ''.join(readl) + '\n')
      if line[0] == '>':
        h = line.strip()
        found_header = True



if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start