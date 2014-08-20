# DO NOT RUN WITH PYTHON2.6! DOESN'T WORK

import sys
import string
import datetime
import random
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import random

def main():
  input_file = sys.argv[1]
  outfile = sys.argv[2]
  plot(input_file, outfile)

def plot(input_file, outfile):
  # Expected format:
  # GGGATAAAGGAGTTT Deg = 13 t = 4 correct

  correct_deg = []
  correct_t = []
  wrong_deg = []
  wrong_t = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    if words[0] != 't-cutoff:':
      if words[7] == 'correct':
        correct_deg.append(int(words[3]))
        correct_t.append(int(words[6]))
      else:
        wrong_deg.append(int(words[3]))
        wrong_t.append(int(words[6]))

  for i in range(len(correct_deg)):
    plt.plot(correct_t[i], correct_deg[i], 'bo', alpha=0.1)

  for i in range(len(wrong_deg)):
    plt.plot(wrong_t[i], wrong_deg[i], 'ro', alpha=0.1)

  plt.xlabel('t of Central Node')
  plt.ylabel('Degree of Central Node')
  # plt.axis([0, 12, 0, 14])
  plt.show()
  # plt.savefig(outfile, format='png')

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start