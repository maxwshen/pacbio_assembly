import sys
import string
import datetime
import random
import copy
from collections import defaultdict

def main():
  if len(sys.argv) != 2:
    print 'Usage: python fastqAnalyze <blasrOutputFile>'
    sys.exit(0)

  fastqAnalyze(sys.argv[1])

def fastqAnalyze(fastq):
  totalbp = 0
  numReads = 0
  isDna = False

  with open(fastq) as f:
    for i, line in enumerate(f):
      if isDna:
        totalbp += len(line)
        numReads += 1
        isDna = False
      if line[0] == '@':
        isDna = True

  print totalbp, numReads

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start