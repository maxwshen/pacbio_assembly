import sys, string, datetime, random, copy, os, fnmatch
import numpy as np
import read_fasta as rf

def main():
  fold = sys.argv[1]
  files = [s for s in os.listdir(fold) if fnmatch.fnmatch(s, '*hood.fasta')]

  for i in range(len(files)):
    print i, len(files)
    fn = files[i]
    reducer(fold, fn)

def reducer(fold, fn):
  step = 5

  with open(fold + fn) as f:
    h, r = rf.read_fasta(fold + fn)

  combined = zip(h, r)
  random.shuffle(combined)
  h[:], r[:] = zip(*combined)

  r = r[:-len(r) % step]   # Make r's len divisible by step
  h = h[:-len(h) % step]
  while len(r) > 0:
    new_f = fold + str(len(r)) + 'x_' + fn
    with open(new_f) as f:
      for i in range(len(r)):
        f.write(h[i] + '\n')
        f.write(r[i] + '\n')
    r = r[:-step]
    h = h[:-step]
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start