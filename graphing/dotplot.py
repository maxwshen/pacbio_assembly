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
  plot(input_file)

def plot(input_file):
  # Expected format:
  # GGGATAAAGGAGTTT Deg = 13 t = 4 correct

  correct_deg = []
  correct_t = []
  wrong_deg = []
  wrong_t = []

  fileh = open(input_file, 'r')
  # for i, line in enumerate(fileh):  
  #   words = line.strip().split()
  #   correctness = words[7]
  #   degree = int(words[3])
  #   _t = int(words[6])
  #   if correctness == 'correct':
  #     correct_deg.append(degree)
  #     correct_t.append(_t)
  #   else:
  #     wrong_deg.append(degree)
  #     wrong_t.append(_t)
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    correctness = words[2]
    degree = int(words[1])
    if correctness == 'correct':
      correct_deg.append(degree)
    else:
      wrong_deg.append(degree)

  largest = [max(correct_deg), max(wrong_deg)]
  binrange = range(0, max(largest))
  for i in range(len(correct_deg)):
    # plt.plot(correct_t[i], correct_deg[i], 'bo', alpha=0.1)
    # plt.plot(0.2, correct_deg[i], 'bo', alpha=0.1)
    pass
  plt.hist(correct_deg, color='green', bins=binrange, alpha=0.4, log=True)

  for i in range(len(wrong_deg)):
    # plt.plot(wrong_t[i], wrong_deg[i], 'ro', alpha=0.1)
    # plt.plot(0, wrong_deg[i], 'ro', alpha=0.1)
    pass
  plt.hist(wrong_deg, color='red', bins=binrange, alpha=0.4, log=True)

  # plt.xlabel('t of Central Node')
  # plt.ylabel('Degree of Central Node')
  plt.xlabel('Degree')
  plt.ylabel('Quantity')
  # plt.axis([0, 30, 0, 10000])
  plt.show()
  # plt.savefig(outfile, format='png')   # doesn't work

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start