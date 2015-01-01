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
  # plot_1group(input_file)
  plot_y_axis(input_file)

def plot_y_axis(input_file):
  fileh = open(input_file, 'r')
  data = fileh.readlines()

  x1 = [i for i in range(len(data))]
  y1 = data

  plt.plot(x1, y1, linewidth = 1.0)
  plt.xlabel('Genome positions')
  plt.ylabel('')
  plt.show()


def plot_1group(input_file):
  x1 = []
  y1 = []

  fileh = open(input_file, 'r')
  data = fileh.readlines()[0].split(', ')
  quantities = []
  quantity = 0
  for i in range(len(data)):
    letter = int(data[i]) + 1
    if len(y1) == 0:
      y1.append(letter)
      quantity += 1
    if len(y1) > 0 and letter == y1[-1]:
      quantity += 1
    if len(y1) > 0 and letter != y1[-1]:
      y1.append(letter)
      quantities.append(quantity)
      quantity = 0

  print y1, '\n', quantities

  x1 = range(1, 34)
  y1 = [3, 15, 4, 5, 14, 5, 22, 5, 18, 12, 9, 5, 19, 2, 21, 20, 3, 15, 13, 13, 5, 14, 20, 19, 19, 15, 13, 5, 20, 9, 13, 5, 19, 4, 15]
  x1 = []
  y1 = []

  for i in range(len(data)):
    x1.append(i)
    y1.append(int(data[i]) + 1)
    print i

  # for i in range(len(x1)):
    # print i
    # plt.plot(x1[i], y1[i], 'go')
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