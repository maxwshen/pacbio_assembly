import sys
import string
import datetime
import random
import copy
from collections import defaultdict

def main():
  if len(sys.argv) != 2:
    print 'Usage: python fastqAnalyze <fastafile>'
    sys.exit(0)

  fastqAnalyze(sys.argv[1])

def fastqAnalyze(fastq):
  totalbp = 0
  numReads = 0
  isDna = False
  counts = dict()

  with open(fastq) as f:
    for i, line in enumerate(f):
      if isDna:
        totalbp += len(line)
        numReads += 1
        isDna = False

        prevChar = ''
        stretch = 0
        for j in range(len(line)):
          if line[j] == prevChar:
            stretch += 1
            if stretch > 3:
              if stretch in counts:
                counts[stretch] += 1
              else:
                counts[stretch] = 1
              for k in range(4, stretch):
                counts[k] += 1
          else:
            stretch = 0
          prevChar = line[j]
      if line[0] == '>':
        isDna = True

  print totalbp, numReads
  print counts

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start