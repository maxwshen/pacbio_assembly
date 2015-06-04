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
  # plot_1group(input_file)
  plot_2groups(input_file)
  # plot_2groups_2files(input_file)

def plot_1group(input_file):
  data = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()

    if float(words[0]) >= 0:
      data.append(float(words[0]))

  largest = max(data)
  # step = int(np.log2(largest)) / 2
  step = 1
  binrange = range(0, largest + step, step)

  # plt.hist(data, color = 'green', bins = binrange, log = True)
  plt.hist(data, color = 'green', bins = binrange)

  print step
  plt.xlabel('Number of (27,6)-mers per read')
  plt.ylabel('Quantity of Reads')
  plt.show()

def plot_2groups(input_file):
  data1 = []
  data2 = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()

    data1.append(float(words[0]))
    data2.append(float(words[1]))

  largest = [max(data1), max(data2)]
  binrange = range(0, max(largest), 3)

  plt.hist(data1, color = 'green', bins = binrange, alpha = 0.4, log = False)
  plt.hist(data2, color = 'red', bins = binrange, alpha = 0.4, log = False)

  plt.xlabel('Initial unc-22 mRNA quantity (ng)')
  plt.ylabel('# Samples')
  plt.show()

def plot_2groups_2files(input_file):
  file1 = 'hist_NEW_acc.txt'
  file2 = 'hist_ORIG_acc.txt'

  data1 = []
  data2 = []

  fileh = open(file1, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    data1.append(float(words[0]))
  
  fileh = open(file2, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()
    data2.append(float(words[0]))    

  largest = [max(data1), max(data2)]
  binrange = range(0, max(largest))

  plt.hist(data1, color = 'green', bins = binrange, alpha = 0.4, log = True)
  plt.hist(data2, color = 'red', bins = binrange, alpha = 0.4, log = True)

  plt.xlabel('Accuracy')
  plt.ylabel('Quantity')
  plt.show()


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start