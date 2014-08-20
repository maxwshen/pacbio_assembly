# Specifically assumes that headers are as follows:
# >read_0/genome_0_10000/mat8523_mis155_ins1068_del254/0_10000
# For example, sim_reads.fasta

import sys
import string
import datetime
import random
import copy
from collections import defaultdict

def main():
  if len(sys.argv) != 5:
    print 'Usage: python section <readsFile> <outFile> <startpos> <endpos>'
    sys.exit(0)

  section(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))

def section(reads, outFile, start, end):
  grabread = False
  with open(outFile, 'w') as out:
    with open(reads) as f:
      for i, line in enumerate(f):
        if grabread:
          grabread = False
          out.write(line)
        if line[0] == '>':
          readstart = int(line.split('/')[1].split('_')[1])
          readend = int(line.split('/')[1].split('_')[2])
          if start <= readstart and readend <= end:
            grabread = True
            out.write(line)


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start