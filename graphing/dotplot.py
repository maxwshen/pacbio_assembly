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
  x1 = []
  y1 = []
  x2 = []
  y2 = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    # correctness = words[2]
    x = int(words[1])
    y = int(words[2])
    # if correctness == 'correct':
    #   x1.append(x)
    #   y1.append(y)
    # else:
    #   x2.append(x)
    #   y2.append(y)
    x1.append(x)
    y1.append(y)

  for i in range(len(x1)):
    plt.plot(x1, y1, 'go', alpha=0.1)

  # for i in range(len(x2)):
    # plt.plot(x2, y2, 'ro', alpha=0.1)

  plt.xlabel('Multiplicity')
  plt.ylabel('Occurences in Genome')
  # plt.axis([0, 30, 0, 10000])
  plt.show()

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start