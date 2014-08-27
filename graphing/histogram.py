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
  data1 = []
  data2 = []

  fileh = open(input_file, 'r')
  for i, line in enumerate(fileh):  
    words = line.strip().split()

    data = int(words[1])
    correctness = words[2]
    
    if correctness == 'correct':
      data1.append(data)
    else:
      data2.append(data)

  largest = [max(correct_deg), max(wrong_deg)]
  binrange = range(0, max(largest))

  plt.hist(data1, color = 'green', bins = binrange, alpha = 0.4, log = True)
  plt.hist(data2, color = 'red', bins = binrange, alpha = 0.4, log = True)

  plt.xlabel('Degree')
  plt.ylabel('Quantity')
  plt.show()

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start