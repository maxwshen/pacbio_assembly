import sys
import string
import datetime
import random
import copy
from collections import defaultdict

def main():
  if len(sys.argv) != 5:
    print 'Usage: python section_sam <input_sam_file> <outFile> <startpos> <endpos>'
    sys.exit(0)

  section_sam(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))

# Filters a sam alignment file by starting/ending position in the genome.
# If an alignment begins within the range, the aligned read is written out
# in .fasta format
def section_sam(inputSam, outFile, startpos, endpos):
  name = ''
  pos = -1
  total = 0
  accepted = 0

  with open(inputSam) as f:
    with open(outFile, 'a') as out:
      for i, line in enumerate(f):
        if line[0] != '@':
          name = str(line.split()[0])
          pos = float(line.split()[3])
          seq = str(line.split()[9])
          total += 1
          if startpos < pos < endpos:
            out.write('>' + name + '\n' + seq + '\n')
            accepted += 1
  print 'Found ' + str(accepted) + ' alignments out of ' + str(total) + '. ' + str(100*float(accepted)/float(total)) + '%'

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start