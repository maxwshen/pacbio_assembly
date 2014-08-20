# findTrueKmer.py
# 
# Takes in reads. Attempts to find "true" kmers by 
# creating a graph where nodes are kmers and edges connect
# kmers that are 1 insertion or 1 deletion away from each
# other.
#
# Implementated with a table approach
#
# Genome file is used to measure accuracy.

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
  genome_sim = '/home/mshen/research/bin/GenomeSimulator.py'

  reads = sys.argv[1]
  genome = sys.argv[2]
  # genome = '/home/mshen/research/sim/sim_genome_1kb_noshift.fasta'
  _k = str( 15 )
  cutoff_deg_nodes = str( 1000 )
  t_greater_than_cutoff = str( 2 )
  _l = str( 0 )
  trim_num = str( 0 )
  num_iterations = 2
  
  t_atleast_cutoff = str(int(t_greater_than_cutoff) + 1)
  # name_base = 'sim_ns1'
  name_base = sys.argv[3]



  for i in range(num_iterations):
    print 'Reads:', reads
    print 'Genome:', genome
    print 'k:', _k
    print 'Cutoff for high degree nodes:', cutoff_deg_nodes
    print 't-cutoff: t >', t_greater_than_cutoff
    print 'l:', _l
    print 'Trim num:', trim_num
    print 'Number of Iterations:', num_iterations, '\n'

    _i = str(i)
    fold = name_base + '_fold_s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '/'
    highdegnodes_outfile = fold + 'highdegnodes.' + _k + '.s' + _i + '.t' + t_atleast_cutoff + '.out'
    highdegnodes_kmers_outfile = fold + 'highdegnodes.' + _k + '.s' + _i + '.t' + t_atleast_cutoff + '.kmers.out'
    contigs_outfile = fold + 'contigsout.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.txt'
    just_contigs_file = fold + 'contigs.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.txt'
    gv_file = fold + name_base + '.s' + _i + '.t' + t_atleast_cutoff + '.L' + _l + '.gv'
    output_reads = name_base + '_substreads.s' + _i + '.t' + t_atleast_cutoff + '.' + _k + '.L' + _l + '.fasta'

    if not os.path.exists(fold):
        os.makedirs(fold)

    # Find high degree nodes in reads
    with open(highdegnodes_outfile, 'w') as f:
      call(python_ver + ' ' + findTrueKmer + ' ' + reads + ' ' + genome + ' ' + _k + ' ' + cutoff_deg_nodes, stdout = f, shell = True)
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

    # Replace contigs into reads
    call(python_ver + ' ' + genome_sim + ' -f ' + fold + ' -c ' + just_contigs_file + ' --kmers ' + highdegnodes_kmers_outfile + ' -g ' + genome + ' -r ' + reads + ' -o ' + output_reads + ' --k ' + _k, shell=True)

    _k = str(int(_k) - 1)
    reads = fold + output_reads

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start