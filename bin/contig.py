# contig.py
# Modified from CSE 181
# 
# Relies on assembly.py and locAL.py
#
# Finds <k> <t> mers in an input fasta file and forms contigs 
# using a debruijn graph assembly approach.
#
# Input file and genome expected to be in fasta format.
# Genome expected to just have 1 sequence

import sys
import string
import datetime
import random
import copy
import assembly
import locAL
import os

from collections import defaultdict

def main():
  if len(sys.argv) < 5:
    print 'Usage: python contig <input file> <genome file> <k> <t> <optional: graphviz name> <optional: l> <optional: trimmingNum>'
    sys.exit(0)

  global kparam
  global tparam
  global lparam
  global trimparam
  kparam = int(sys.argv[3])
  tparam = int(sys.argv[4])
  gvname = None
  lparam = None
  trimparam = 0
  try:
    trimparam = int(sys.argv[7])
    gvname = sys.argv[5]
    lparam = int(sys.argv[6])
    if lparam >= kparam:
      print 'ERROR: Bad l: Greater than or equal to k'
      sys.exit()
    if trimparam < 1:
      trimparam = 0
  except:
    pass

  # kmers = assembly.findKTmers(sys.argv[1], kparam, tparam)
  # if lparam is not None:
  #   kmers = splitKTmers(kmers, kparam - lparam)

  # Read kmers in manually from high degree nodes
  kmers = set()
  with open('findTrueKmers.nhood_sim4.10885.15.top1k.txt') as f:
    for i, line in enumerate(f):
      kmers.add(line.strip())
  # print kmers

  numEdges = 0
  adjlist, numEdges = deBruijnkmers(kmers)

  adjlist = trimdebruijn(adjlist, trimparam)
  contig(adjlist, sys.argv[2], gvname)


# Input:
#   Adjacency list of a de bruijn graph
# Output:
#   Trimmed adjacency list of a de bruijn graph
# Procedure:
#   Remove all nodes with zero indegree, iterate as desired
def trimdebruijn(adjlist, num):
  for i in range(num):
    degrees = dict()
    outdegree = dict()
    indegree = dict()
    inedges = []
    for node in adjlist:
      degrees[node] = len(adjlist[node])
      outdegree[node] = len(adjlist[node])
      indegree[node] = 0
    for suffixlists in adjlist.values():
      inedges += suffixlists
    for suffix in inedges:
      if suffix in degrees:
        degrees[suffix] -= 1
      else:
        degrees[suffix] = -1
      if suffix in indegree:
        indegree[suffix] += 1
      else:
        indegree[suffix] = 1
      if suffix not in outdegree:
        outdegree[suffix] = 0

    trimcount = 0
    # Get node w/ 0 indegree
    for k, v in indegree.iteritems():
      if v == 0:
        print ' ', k, indegree[k], outdegree[k], adjlist[k]
        del adjlist[k]
        trimcount += 1

    # Get node w/ 0 outdegree
    for k, v in outdegree.iteritems():
      if v == 0:
        for suffixlist in adjlist.values():
          if k in suffixlist:
            suffixlist.remove(k)
            trimcount += 1
            break

    print 'Trimmed', trimcount, 'leaves. New graph size:', len(adjlist), 'nodes'
  return adjlist

def contig(adjlist, genomeFile, gvname):
  if gvname == None:
    gvname = 'graph_viz.gv'
  # for i in adjlist:
  #   print i, adjlist[i]
  # print '---'
  # print numEdges

  # Construct degree list for all nodes to find where to start. 
  # Should start at node that has total outdegree of 1
  degrees = dict()
  outdegree = dict()
  indegree = dict()
  inedges = []
  for node in adjlist:
    degrees[node] = len(adjlist[node])
    outdegree[node] = len(adjlist[node])
    indegree[node] = 0
  for suffixlists in adjlist.values():
    inedges += suffixlists
  for suffix in inedges:
    if suffix in degrees:
      degrees[suffix] -= 1
    else:
      degrees[suffix] = -1
    if suffix in indegree:
      indegree[suffix] += 1
    else:
      indegree[suffix] = 1
    if suffix not in outdegree:
      outdegree[suffix] = 0
  
  # print degrees
  # print outdegree
  # print indegree
  # print '----'

  # Get node w/ total outdegree of 1
  for k, v in degrees.iteritems():
    if v == 1:
      # print ' ', k
      begin = k

  cycle = []
  nodes = adjlist.keys() # prefixes
  cycle.append(begin)
  current = begin
  intermediate = []
  contigs = []
  noutdegree = copy.deepcopy(outdegree)
  nindegree = copy.deepcopy(indegree)

  try:
    open(gvname, 'w').close()
  except:
    pass
  with open(gvname, 'a') as f:
    f.write('digraph G {\n')
    while True:
      #print current, adjlist[current][0], contigs
      f.write('\t' + current + ' -> ' + adjlist[current][0] + ';\n')

      intermediate.append(current)
      prev = current
      if len(adjlist[current]) == 0:
        # print 'breaking on', current, adjlist[current]
        break
      current = adjlist[current][0]
      adjlist[prev] = adjlist[prev][1:]   # Remove edge we just traversed
      noutdegree[prev] -= 1
      nindegree[current] -= 1
      # print '  out', outdegree
      # print '  in', indegree
      if outdegree[current] == 1 and indegree[current] == 1:
        continue
      else: # Found branching point. Put intermediate contig into contigs
        # print 'stopped'
        intermediate.append(current)
        contigs.append(intermediate)
        intermediate = []
        # If we're at the end of a path, start from another node w/ 0 indegree
        if len(adjlist[current]) == 0:
          # print 'TEST'
          current = -1
          for n in nodes:
            if nindegree[n] == 0 and noutdegree[n] > 0:
              current = n
          if current == -1:
            break
    f.write('}')

  # Ensure that graph is fully traversed
  for i in adjlist:
    if len(adjlist[i]) > 0:
      print 'ERROR:', i, 'still has untraversed edges'

  # Put contigs together
  contigs.sort()
  c = []
  for contig in contigs:
    temp = contig[0]
    for i in contig[1:]:
      temp += i[-1:]
    c.append(temp)

  totalbp = 0
  for contig in c:
    print contig
    totalbp += len(contig)

  # return

  # Check the accuracy of the contigs by aligning to the genome file
  print 'Aligning to genome...'
  with open(genomeFile) as f:
    lines = f.readlines()
    genome = lines[1]

  numPerfect = 0
  matchScores = []
  perfectStarts = []
  perfectLens = []
  for i in range(len(c)):
    print c[i]
    (alignLen, matches, mismatches, numgaps, numGapExtends, bestxy) = locAL.external_bestseq1(c[i], genome, 1, -1, -1, -0.5)
    print 'Range:', bestxy[1] - alignLen, '-', bestxy[1]
    if alignLen != 0:
      score = float(matches)/float(alignLen)
      print 'ERROR: Divide by Zero! Alignment Length = 0'
    else:
      score = 0
    matchScores.append(score)
    if score == 1:
      numPerfect += 1
      perfectStarts.append(bestxy[1] - alignLen)
      perfectLens.append(alignLen)
  perfectLens = [x for (y,x) in sorted(zip(perfectStarts, perfectLens))]
  perfectStarts = sorted(perfectStarts)

  # for score in matchScores:
  #   print str(score) + '\t',
  # print '\n'

  perfectRanges = [] # stores tuples
  start = -1
  end = -1
  for i in range(len(genome)):
    if i in perfectStarts:
      if end < i + perfectLens[perfectStarts.index(i)]:
        end = i + perfectLens[perfectStarts.index(i)]
      if start == -1:
        start = i
    if i == end:
      perfectRanges.append((start, end))
      start = -1
      end = -1

  totalCovered = 0
  for (x, y) in perfectRanges:
    totalCovered += y - x

  print 'k =', kparam, ', t =', tparam, ', l =', lparam, ', #trims =', trimparam
  print '# contigs:', len(c)
  print '# bp total:', totalbp
  print '# avg contig len:', float(totalbp)/float(len(c))
  print '# perfect:', numPerfect, '\t', float(numPerfect*100)/float(len(c)), '%'
  print 'Avg. accuracy:', 100*sum(matchScores)/len(matchScores), '%'
  print 'Perfectly Covered Ranges:', perfectRanges
  print 'Total covered perfectly:', totalCovered, '\tPercent:', float(totalCovered)/float(len(genome))


# Input:
#   A list of kmers
# Output:
#   deBruijn adjacency list for each kmer using prefix and suffix,
#   formatted as a dictionary of prefix keys and suffix values.
#   If there are multiple edges between the same prefix and suffix,
#   there will be multiple suffixes in the dictionary. Example:
#   dict['ACTG'] = ['CTGG', 'CTGG']
# Process:
#   Separate all kmers into prefixes and suffixes of length (k-1)
#   Make an edge between prefixes and suffixes that overlap
def deBruijnkmers(kmers_in):
  numEdges = 0
  prefixes = []
  suffixes = []
  for kmer in kmers_in:
    prefixes.append(kmer[:len(kmer)-1])
    suffixes.append(kmer[1:])

  kmers = defaultdict(list)
  for k in prefixes:
    kmers[k].append('')

  for i in range(len(prefixes)):
    prefix = prefixes[i]
    suffix = suffixes[i]
    if overlaps(prefix, suffix):
      numEdges += 1
      kmers[prefix].append(suffix)

  for k in prefixes:
    kmers[k].remove('')

  return kmers, numEdges

# Returns true if suffix of q1 overlaps with the prefix of q2
# Used in deBruijnkmers
def overlaps(q1, q2):
  q1 = q1[1:]
  q2 = q2[:len(q2)-1]
  if len(q1) != len(q2):
    return False
  for i in range(len(q1)):
    if q1[i] != q2[i]:
      return False
  return True

# Input:
#   Set of kmers
# Output:
#   Set of l-mers created from the kmers 
def splitKTmers(kmers, l):
  print 'Splitting ktmers with l =', l
  newkmers = set()
  for kmer in kmers:
    for i in range(len(kmer) - l + 1):
      newkmers.add(kmer[i:i+l])
  return newkmers


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start