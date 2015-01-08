import sys, string, datetime, random, copy, os, commands, fnmatch, numpy
from collections import defaultdict

import read_fasta as rf

def main():
  # genomic_coverage_abruijn_reads()
  infile = sys.argv[1]
  get_coverage(infile)
  return

def get_coverage(infile):
  # inp = '/home/mshen/research/nohup_gen_coverage_abruijn_reads_22.4.out'
  inp = infile
  with open(inp) as f:
    aligns = f.readlines()
  genome = [0] * 4700000
  for a in aligns:
    if len(a.split()) > 7:
      beg = int(a.split()[6])
      end = int(a.split()[7])
      for i in range(beg, end):
        genome[i] += 1

  regions_0 = []
  start_pos = -1
  extend = False
  for i in range(len(genome)):
    if genome[i] == 0 and not extend:
      extend = True
      start_pos = i
    if genome[i] != 0 and extend:
      extend = False
      regions_0.append([start_pos, i])

  # print regions_0
  lens_0 = [s[1] - s[0] for s in regions_0]
  print lens_0
  print numpy.mean(lens_0), numpy.std(lens_0), len(regions_0)
  print sorted(regions_0, key = lambda d: d[1] - d[0], reverse = True)
  sys.exit(0)

  # out_file = 'out_abruijn_coverage_22.4.out'
  out_file = 'out_abruijn_contigs4.out'
  with open(out_file, 'w') as f:
    f.write('\n'.join([str(s) for s in genome]))
  return

def genomic_coverage_abruijn_reads():
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'   # Concise output
  temp_sig = str(datetime.datetime.now()).split()[1]

  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)

  traversed_headers = set()
  temp_file = 'temp_' + temp_sig + '.fasta'
  for kt in headers.keys():
    for h in headers[kt]:
      if h not in traversed_headers:
        traversed_headers.add(h)
        with open(temp_file, 'w') as f:
          f.write(h + '\n' + rr[hr.index(h)])
        status = commands.getstatusoutput(blasr_exe + ' ' + temp_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]
        if len(status) != 0:
          print status

def build_headers_dict(ktmer_headers_file):
  headers = defaultdict(list)   # Key = ktmer, Val = [headers]
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      headers[words[0]] = words[1:]
  return headers

if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start