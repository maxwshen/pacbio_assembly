# Run with python2.4 on debruijn server

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
  plot_1group(input_file)

def plot_1group(input_file):
  x1 = []
  y1 = []

  fileh = open(input_file, 'r')
  data = fileh.readlines()[0].split(', ')
  for i in range(len(data)):
    x1.append(i)
    y1.append(int(data[i]))
    print i

  for i in range(len(x1)):
    print i
    plt.plot(x1, y1, 'go')

  plt.xlabel('time')
  plt.ylabel('likeliest hidden state')
  plt.show()

def plot_2groups(input_file):
  x1 = []
  y1 = []
  x2 = []
  y2 = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    correctness = words[7]
    x = int(words[6])
    y = int(words[3])
    if correctness == 'correct':
      x1.append(x)
      y1.append(y)
    else:
      x2.append(x)
      y2.append(y)
    # x1.append(x)
    # y1.append(y)

  for i in range(len(x1)):
    plt.plot(x1, y1, 'go', alpha=0.1)

  for i in range(len(x2)):
    plt.plot(x2, y2, 'ro', alpha=0.1)

  plt.xlabel('Multiplicity')
  plt.ylabel('Degree')
  # plt.axis([0, 30, 0, 10000])
  plt.show()

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start