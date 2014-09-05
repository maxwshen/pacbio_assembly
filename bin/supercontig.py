import sys
import string
import datetime
import random
import copy
import os
import numpy as np

from collections import defaultdict
from subprocess import call

def main():
  python_ver = 'python2.6'
  findTrueKmer = '/home/mshen/research/bin/findTrueKmer.py'
  contig_nodes = '/home/mshen/research/bin/contig_nodes.py'
  assemble_contigs = '/home/mshen/research/bin/assemble_contigs.py'
  genome_sim = '/home/mshen/research/bin/GenomeSimulator.py'

  reads = sys.argv[1]
  genome = sys.argv[2]
  # genome = '/home/mshen/research/sim/sim_genome_1kb_noshift.fasta'
  _k = str( 15 )
  _d = str( 2 )
  cutoff_deg_nodes = str( 1000 )
  t_greater_than_cutoff = str( 1 )
  filter_neighbors = 'True'
  _l = str( 0 )
  trim_num = str( 0 )
  
  t_atleast_cutoff = str(int(t_greater_than_cutoff) + 1)
  # name_base = 'sim_ns1'
  name_base = sys.argv[3]


  for j in range(10, 9, -1):
    for i in range(0, 1):
      _i = str(i)
      _j = str(j)

      # Simulated reads
      # j range from 30 to 5, i range 10
      # reads = '/home/mshen/research/sim/sim_10ns/sim_reads_1kb_ns' + _i + '_cov' + _j + '.fasta'
      # genome = '/home/mshen/research/sim/sim_10ns/sim_genome_1kb_ns' + _i + '.fasta'

      # E. coli reads
      # j range from 35 to 5, i range 12
      # reads = '/home/mshen/research/data/all_cov/ec_reads_rh_hc_n' + _i + '_cov' + _j + '.fasta'
      # genome = '/home/mshen/research/data/all_cov/ec_genome_rh_hc_n' + _i + '.fasta'

      # print 'Reads:', reads
      # print 'Genome:', genome
      # print 'k:', _k
      # print 'Cutoff for high degree nodes:', cutoff_deg_nodes
      # print 't-cutoff: t >', t_greater_than_cutoff
      # print 'l:', _l
      # print 'Trim num:', trim_num, '\n'

      # fold = name_base + '_fold_s' + _i + '.cov' + _j + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '/'
      fold = name_base + '_fold_s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '/'
      highdegnodes_outfile = fold + 'highdegnodes.' + _k + '.s' + _i + '.t' + t_atleast_cutoff + '.out'
      highdegnodes_kmers_outfile = fold + 'highdegnodes.' + _k + '.s' + _i + '.t' + t_atleast_cutoff + '.kmers.out'
      contigs_outfile = fold + 'contigsout.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.txt'
      just_contigs_file = fold + 'contigs.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.txt'
      gv_file = fold + name_base + '.s' + _i + '.t' + t_atleast_cutoff + '.L' + _l + '.gv'
      output_reads = name_base + '_substreads.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.fasta'
      assemble_outfile = fold + 'assembly.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.txt'

      if not os.path.exists(fold):
          os.makedirs(fold)

      # Find high degree nodes in reads
      with open(highdegnodes_outfile, 'w') as f:
        call(python_ver + ' ' + findTrueKmer + ' ' + reads + ' ' + genome + ' ' + _k + ' ' + _d + ' ' + cutoff_deg_nodes + ' ' + t_atleast_cutoff + ' ' + filter_neighbors, stdout = f, shell = True)
      with open(highdegnodes_outfile) as f:
        with open(highdegnodes_kmers_outfile, 'w') as g:
          for i, line in enumerate(f):
            if line.split()[0] != 't-cutoff:':
              g.write(line)
            else:
              print line

      # Find super-contigs from high degree kmers by building de Bruijn graph
      with open(contigs_outfile, 'w') as f:
        call(python_ver + ' ' + contig_nodes + ' ' + highdegnodes_kmers_outfile + ' ' + genome + ' ' + _k + ' ' + t_greater_than_cutoff + ' ' + gv_file + ' ' + _l + ' ' + trim_num, stdout = f, shell = True)

      get_contigs = False
      with open(contigs_outfile) as f:
        with open(just_contigs_file, 'w') as g:
          for i, line in enumerate(f):
            if line == '\n':
              get_contigs = False
            if get_contigs:
              g.write(line)
            if line == 'Done\n':
              get_contigs = True

      # Assemble super-contigs
      with open(assemble_outfile, 'w') as f:
        call(python_ver + ' ' + assemble_contigs + ' ' + contigs_outfile + ' ' + highdegnodes_kmers_outfile + ' ' + genome, stdout = f, shell =  True)

      print reads
      with open(assemble_outfile) as f:
        for i, line in enumerate(f):
          if line.split()[0] == 'Range:' or line.split()[0] == 'Accuracy:':
            print line


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start