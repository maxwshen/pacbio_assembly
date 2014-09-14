import sys
import string
import datetime
import random
import copy
from collections import defaultdict

def main():
  blasrOutputFile = sys.argv[1]
  readsFile = sys.argv[2]
  outFile = sys.argv[3]
  startpos = int(sys.argv[4])
  endpos = int(sys.argv[5])
  threshold = float(sys.argv[6])

  section(blasrOutputFile, readsFile, outFile, startpos, endpos, threshold)
  return

def section(blasrOut, reads, outFile, start, end, thresholdPct):
  threshold = thresholdPct
  valid = set()

  # Put query names for reads aligned in range and >50% aligned
  # into the set 'valid' 
  with open(blasrOut) as f:
    query = ''
    tstart = -1
    tend = -1
    for i, line in enumerate(f):
      if len(line) > 10:
        # If it's the query name
        if line[9] == 'Q':
          query = str(line.split()[1])
          query = query[:query.rfind('/')]
        # If it's the query range
        if line[3] == 'Q':
          ln = line.split()
          qstart = float(ln[1])
          qend = float(ln[3])
          qtotal = float(ln[5])
          qpct = (qend - qstart) / qtotal
        # If it's the target range
        if line[2] == 'T':
          ln = line.split()
          tstart = float(ln[1])
          tend = float(ln[3])
        if start < tstart < tend < end and qpct > threshold:
          valid.add(query)
          query = ''
          qpct = 0
          tstart = -1
          tend = -1

  # Using the set of valid query names, filter reads
  found = False
  with open(reads) as f:
    for i, line in enumerate(f):
      if line[0] == '@':
        found = False
        if str(line[1:].strip()) in valid:
          found = True
      if found:
        with open(outFile, 'a') as fout:
          fout.write(line)

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start