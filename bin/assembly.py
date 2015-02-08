# A-Bruijn Graph Construction

import sys, string, datetime, random, copy, os
import numpy as np
from collections import defaultdict
import find_read

def main():
  # reads_file = sys.argv[1]
  # genome_file = sys.argv[2]
  # _k = int(sys.argv[3])
  # _t = int(sys.argv[4])
  # gv_file = sys.argv[5]
  reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  # reads_file = '/home/mshen/research/sample.fasta'
  genome_file = '/home/mshen/research/data/e_coli_genome.fasta'
  _k = int(sys.argv[1])
  _t = int(sys.argv[2])
  gv_file = 'temp.gv'
  assembly(reads_file, genome_file, _k, _t, gv_file)
  return

def dist_bw_ktmers_in_genome(ktmers, _k, genome_file, out_file):
  genome = ''
  with open(genome_file) as f:
    for i, line in enumerate(f):
      if i > 0:
        genome += line.strip()

  dists = []
  curr = 0
  for i in range(len(genome) - _k + 1):
    curr += 1
    if genome[i : i + _k] in ktmers:
      dists.append(curr)
      curr = 0

  dists.append(curr)
  with open(out_file, 'w') as f:
    f.write('\n'.join([str(s) for s in dists]))
  return


def assembly(reads, genome_file, _k, _t, gvname):
  # Given input reads in fasta format, assemble the genome
  # Procedure:
  #   1. Find all (k,T)-mers in reads
  #   2. Represent reads as list of alternating (k,T)-mers
  #   3. Construct A-Bruijn graph with reads

  print 'Assembling', reads, 'with k =', _k, 'and t =', _t
  print 'Finding kt-mers...', 
  ktmers = findKTmers(reads, _k, _t)
  print '... Done.', datetime.datetime.now()

  # Find distance bw kt-mers for figures
  out_file = 'dists_' + str(_k) + '_' + str(_t) + '_genome.txt'
  dist_bw_ktmers_in_genome(ktmers, _k, genome_file, out_file)

  print 'Converting reads...'
  cReads = convertReads(reads, ktmers, _k)
  print '... Done.', datetime.datetime.now()

  headers_out_file = 'temp_ktmer_headers' + '_' + str(_k) + '_' + str(_t) + '.out'
  edges_out_file = 'temp_ktmer_edges' + '_' + str(_k) + '_' + str(_t) + '.out'
  creads_out_file = 'temp_creads.out' + '_' + str(_k) + '_' + str(_t) + '.out'
  a_bruijn_summary(cReads, reads, headers_out_file, edges_out_file, creads_out_file)
  return

  # # Writes distance statistics to file
  # distfn = 'diststats_' + reads + '.' + str(_k) + '.' + str(_t) + '.txt'
  # distfn = distfn.translate(None, '/')
  # try:
  #   os.remove(distfn)
  # except:
  #   pass
  # with open(distfn, 'a') as f:
  #   for i in cReads:
  #     for item in i[1]:
  #       if item > 0:
  #         f.write(str(item) + '\t')

  print 'Generating graph...', datetime.datetime.now()
  graph = A_Bruijn_Graph(cReads)
  print 'Total number of nodes:', len(graph.allnodes)
  print '... Done.', datetime.datetime.now()

  print 'Splitting large nodes...'
  graph.split_large_nodes(reads)
  print 'Total number of nodes after splitting:', len(graph.allnodes)
  print '... Done.', datetime.datetime.now()

  print 'Finding neighborhoods...'
  genome = ''
  with open(genome_file) as f:
    for i, line in enumerate(f):
      if i > 0:
        genome += line.strip()

  starting_pos = 250          # 625
  jump_length = 400           # 1250
  neighborhood_width = 500    # 1250
  neighborhood_margin = 500   # 0
  position = starting_pos
  while position < len(genome) - _k + 1:
    seedktmer = genome[position : position + _k]
    if seedktmer in [g[:_k] for g in graph.allnodes]:
      for g in graph.allnodes:
        if genome[position : position + len(g)] == g:
          curr_ktmer = g
      seednode = graph.allnodes[curr_ktmer]
      neighborhood(reads, seednode, neighborhood_width, neighborhood_margin, position, graph)
      print '... Done.', position, datetime.datetime.now()
      position += jump_length
    else:
      position += 1

  print '... Done.', datetime.datetime.now()

  return


  print 'Traversing graph...'
  bulges = 0
  samektmercount = 0
  num1edge = 0
  num1edge1 = 0
  num1edge1in1 = 0
  totalkmerlen = 0
  numCombinedNodes = 0
  totalCombkmerlen = 0
  traversed = set()

  nodes = set()
  for n in headNode.outnodes:
    nodes.add(n)
  
  try:
    open(gvname, 'w').close()
  except:
    pass
  with open(gvname, 'a') as f:
    f.write('digraph G {\n')
    while len(nodes) > 0:
      # print len(nodes)
      n = nodes.pop()
      traversed.add(n)
      totalkmerlen += len(n.ktmer)
      if len(n.ktmer) > _k:
        totalCombkmerlen += len(n.ktmer)
        numCombinedNodes += 1
      # print n.ktmer, n.outDegree, n.inDegree
      if n.inDegree > 0: 
        bulges += n.inDegree - 1
      print n.ktmer + ':', n.outDegree, 'edges out and', n.inDegree, 'edges in'
      findSameKtmers = set()
      for i in range(len(n.outedges)):

        # Graphing
        f.write('\t' + n.ktmer + ' -> ' + n.outnodes[i].ktmer + ' [label = "' + str(n.outedges[i]) + '"];\n')

        if n.outnodes[i].ktmer in findSameKtmers:
          samektmercount += 1
        else:
          findSameKtmers.add(n.outnodes[i].ktmer)
        print '  dist:', n.outedges[i], n.outnodes[i].ktmer, '\treadcount:', n.outReadCount[i]
        if n.outnodes[i] not in traversed:
          nodes.add(n.outnodes[i])
      if n.outDegree == 1:
        num1edge += 1
        if n.outedges[0] == 1:
          num1edge1 += 1
        if n.outnodes[0].inDegree == 1:
          num1edge1in1 += 1
    f.write('}')

  print 'Assembled', reads, 'with k =', _k, 'and t =', _t
  print '# Nodes joined:', numJoined
  print '# Bulges:', bulges
  print '# Same edge, different dist:', samektmercount 
  print '# Nodes:', len(traversed)
  print '# Nodes w/ exactly 1 out edge:\t\t', num1edge, str(float(num1edge*100)/float(len(traversed))) + '%'
  print '# Nodes w/ exactly 1 out edge, dist 1:\t', num1edge1, str(float(num1edge1*100)/float(len(traversed))) + '%'
  print '# Edges b/w nodes of outdegree/indegree 1:\t', num1edge1in1, str(float(num1edge1in1*100)/float(len(traversed))) + '%'
  print 'Avg. kt-mer len:', str(float(totalkmerlen)/float(len(traversed)))
  
  if numCombinedNodes > 0:
    print '# Combined Nodes:', numCombinedNodes
    print 'Avg. kt-mer len for combined nodes:', str(float(totalCombkmerlen)/float(numCombinedNodes))
  else:
    print '# Combined Nodes:', numCombinedNodes
    print 'Avg. kt-mer len for combined nodes: N/A'


def a_bruijn_summary(cReads, reads_file, headers_out_file, edges_out_file, creads_out_file):
  # Make a defaultdict(list) of edges for each node, and distance
  with open(creads_out_file, 'w') as f:
    for c in cReads:
      line = [c[1][0]]
      for i in range(len(c[0])):
        line.append(c[0][i])
        line.append(c[1][i + 1])
      f.write(' '.join([str(s) for s in line]) + '\n')

  headers = dict()  # Key = number, Value = header
  num = 0
  with open(reads_file) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        headers[num] = line.strip()
        num += 1

  reads_kt = defaultdict(list)  # Key = ktmer, Value = [header1, dist_to_start, dist_to_end, header2, ...]
  edges = defaultdict(list)     # Key = ktmer, Value = [(neighbor, dist), ...]

  minimum = 20
  for i in range(len(cReads)):
    if i % 5000 == 0:
      print i, datetime.datetime.now()
    # print cReads[i]
    h = headers[i]
    ktmers = cReads[i][0]
    dists = cReads[i][1]
    for j in range(len(ktmers)):
      kt = ktmers[j]
      if h not in reads_kt[kt]:
        reads_kt[kt].append(h.split()[0])
        # r = find_read.find_read(h, reads_file)
        # r = r.splitlines()[1].strip()
        # reads_kt[kt].append(r.index(kt))
        # reads_kt[kt].append(len(r) - r.index(kt) - 1)
      # for k in range(len(ktmers)):
      #   if k < j and ktmers[k] not in [s[0] for s in edges[kt]]:
      #     dist = -1 * sum(dists[k + 1 : j + 1])
      #     if abs(dist) > minimum:
      #       edges[kt].append((ktmers[k], dist))
      #       edges[ktmers[k]].append((kt, - dist))
      #   if k > j and ktmers[k] not in [s[0] for s in edges[kt]]:
      #     dist = sum(dists[j + 1 : k + 1])
      #     if abs(dist) > minimum:
      #       edges[kt].append((ktmers[k], dist))
      #       edges[ktmers[k]].append((kt, - dist))

  with open(headers_out_file, 'w+') as f:
    for k in reads_kt.keys():
      f.write(k + ' ' + ' '.join([str(s) for s in reads_kt[k]]) + '\n')

  # with open(edges_out_file, 'w+') as f:
  #   for k in edges.keys():
  #     out_edges = [s[0] + ' ' + str(s[1]) for s in edges[k]]
  #     f.write(k + ' ' + ' '.join(out_edges) + '\n')

  return

def neighborhood(reads, centerNode, dist, margin, position, graph):
  # Input:
  #   node: A node in the A-bruijn graph
  #   dist: Integer
  # Output:
  #   A set of sequences corresponding to the neighborhood of dist/2 around node
  # Margin refers to searching for nodes within the neighborhood by
  #   exploring farther outside and then back in.
  # Position is the true genomic position.
  # allnodes is a dictionary of all nodes, Keys = kmers, Values = Nodes

  allnodes = graph.allnodes
  nodes = dict()        # Keys are ktmers, values are position, 0 = starting pt
  traversed = set()     # Set of ktmers
  collected = dict()    # Keys are ktmers, values are position, 0 = starting pt
  distmin = (dist * -1)/2
  distmax = dist / 2

  current = centerNode
  currpos = 0

  nhood_size_limit = 100
  while True and len(collected) < nhood_size_limit:
    # Add all neighboring nodes if not traversed already
    for i in range(len(current.outnodes)):
      nextNode = current.outnodes[i]
      if nextNode.ktmer not in traversed and nextNode not in nodes.keys():
        nodes[nextNode.ktmer] = currpos + current.outedges[i]
    for i in range(len(current.innodes)):
      nextNode = current.innodes[i]
      if nextNode.ktmer not in traversed and nextNode not in nodes.keys():
        nodes[nextNode.ktmer] = currpos - current.inedges[i]
    
    # Remove nodes that are outside of alloted distance
    tempdict = dict()
    for key in nodes:
      if distmin - margin < nodes[key] < distmax + margin:
        tempdict[key] = nodes[key]
    nodes = tempdict

    traversed.add(current.ktmer)
    if distmin < currpos < distmax:
      collected[current.ktmer] = currpos

    if len(nodes) == 0:
      break
    else:
      curr_ktmer = nodes.keys()[0]
      current = allnodes[curr_ktmer]
      currpos = nodes[curr_ktmer]
      del nodes[curr_ktmer]

    # print len(nodes), '\tNodes:', nodes, '\tCollected:',collected
    print len(nodes)


  # Find starting position and read for regions in neighborhood
  hoodReads = []
  hoodStartPos = []
  for n_tag in collected:
    n = allnodes[n_tag]
    print n.ktmer, n.reads, n.pos, collected[n.ktmer]
    for i in range(len(n.reads)):
      if n.reads[i] not in hoodReads:
        hoodReads.append(n.reads[i])
        hoodStartPos.append(n.pos[i] + distmin - collected[n.ktmer])

  print hoodReads, hoodStartPos

  readNum = -1
  seqs = []
  with open(reads) as f:
    for i, line in enumerate(f):
      if line[0] == '>':
        readNum += 1
        if readNum in hoodReads:
          seqs.append(line.strip())
      if readNum in hoodReads and line[0] != '>':
        startpos = hoodStartPos[hoodReads.index(readNum)]
        # print readNum, line
        if startpos > 0:
          seqs.append(line[startpos:startpos + dist].strip())
        else:
          seqs.append(line[0:startpos + dist].strip())

  fold = '/home/mshen/research/nhoods_split_24.4/'
  if not os.path.exists(fold):
    os.makedirs(fold)
  nhood_filename = fold + 'nhood_nh' + str(position) + '_' + centerNode.ktmer + '.fasta'
  try:
    open(nhood_filename, 'w').close()
  except:
    pass
  with open(nhood_filename, 'a') as f:
    for seq in seqs:
      f.write(seq + '\n')


def findKTmers(reads, _k, _t):
  # Input:
  #   file: Reads in fasta format
  #   int: _k
  #   int: _t
  # Output:
  #   Set of all ktmers in reads  

  counts = dict()
  readcount = 0
  with open(reads) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] != '>':
        dna = line.strip()
        curr_read += dna
      if line[0] == '>' or line[0] == '@':
        readcount += 1
        for j in range(len(curr_read) - _k + 1):
          kmer = curr_read[j:j+_k]
          # print kmer
          if kmer in counts:
            counts[kmer] = counts[kmer] + 1
          else:
            counts[kmer] = 1
        curr_read = ''
    for j in range(len(curr_read) - _k + 1):
      kmer = curr_read[j:j+_k]
      # print kmer
      if kmer in counts:
        counts[kmer] = counts[kmer] + 1
      else:
        counts[kmer] = 1
    curr_read = ''
  ans = set()

  for key, val in counts.iteritems():
    if _t <= val:
      if not_palindrome(key):
        ans.add(key)
  print 'Found', len(ans), 'ktmers in', readcount, 'reads'
  return ans


def not_palindrome(kmer):
  if len(kmer) % 2 == 1:
    if kmer[: len(kmer) / 2 + 1][::-1] == kmer[len(kmer) / 2 :]:
      return False
    return True 
  else:
    if kmer[: len(kmer) / 2][::-1] == kmer[len(kmer) / 2 :]:
      return False
    return True


def convertReads(reads, ktmers, _k):
  # Input:
  #   file: Reads in fasta format
  #   set: ktmers
  #   int: _k
  # Output:
  #   A list of tuples, each tuple is 2 lists. Each tuple describes 1 read.
  #   Each tuple contains two lists, the first list is all (k,t)-mers
  #   and the second list is the distance b/w each (k,t)-mer.
  #   The i-th element in the distance list describes the distance b/w
  #     the (i-1)th (k,t)-mer and the i-th (k,t)-mer.
  #   The 0-th element in the distance list is the distance of the first (k,t)-mer
  #     from the beginning of the string.

  cReads = []
  readcount = 0
  num_reads_wo_ktmers = 0
  curr_dna = ''

  readcount = 0
  with open(reads) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] != '>':
        dna = line.strip()
        curr_read += dna
      if i > 0:
        if line[0] == '>' or line[0] == '@':
          readcount += 1
          tempDist = []
          tempKmers = []
          count = 0
          for j in range(len(curr_read) - _k + 1):
            kmer = curr_read[j : j + _k]
            if kmer in ktmers:
              tempDist.append(count)
              tempKmers.append(kmer)
              count = 1
            else:
              count += 1
          tempDist.append(count + _k - 1)
          if len(tempDist) == 1:
            num_reads_wo_ktmers += 1
          cReads.append((tempKmers, tempDist))
          curr_read = ''
    readcount += 1
    tempDist = []
    tempKmers = []
    count = 0
    for j in range(len(curr_read) - _k + 1):
      kmer = curr_read[j : j + _k]
      if kmer in ktmers:
        tempDist.append(count)
        tempKmers.append(kmer)
        count = 1
      else:
        count += 1
    tempDist.append(count + _k - 1)
    if len(tempDist) == 1:
      num_reads_wo_ktmers += 1
    cReads.append((tempKmers, tempDist))
    curr_read = ''

  # print cReads
  print num_reads_wo_ktmers, 'reads without any kt-mers out of', readcount, 'reads'
  print cReads
  return cReads

def hamming_dist(s1, s2):
  if len(s1) != len(s2):
    return -1
  return sum([s1[i] != s2[i] for i in range(len(s1))])

def kmer_matching(ec_seq, hr, rr, _k, cutoff):
  if len(hr) == 0 or len(rr) == 0:
    return None
  kmers = set()
  # for i in range(len(ec_seq) - _k + 1):
  for i in range((len(ec_seq) - 500) / 2, (len(ec_seq) - 500) / 2 + 500):
    kmers.add(ec_seq[i:i + _k])

  reads = dict()    # Key = header, value = num shared kmers
  for i in range(len(rr)):
    r = rr[i]
    h = hr[i]
    score = sum([1 if r[i:i + _k] in kmers else 0 for i in range(len(r) - _k + 1)])
    reads[h] = score

  headers = []
  for key in sorted(reads, key = reads.get, reverse = True):
    if reads[key] < cutoff:
      break
    headers.append(key)
    print 'km score:', reads[key]
  return headers

def split_2_groups(reads):
  # Splits reads from kmer matching into 2 groups based on minimizing stdev
  # Doesn't exactly work well with kmer matching as it is - stdev not a good divider
  best_score = float('inf')
  bestg1 = []
  for i in range(1, len(reads)):
    group1 = []
    group2 = []
    curr_pos = 0
    for key in sorted(reads, key = reads.get, reverse = True):
      if curr_pos < i:
        group1.append(key)
      else:
        group2.append(key)
      curr_pos += 1
    score = np.std([reads[key] for key in group1]) + np.std([reads[key] for key in group2])
    if score < best_score:
      bestg1 = group1
      best_score = score
  print best_score, bestg1
  for k in reads.keys():
    print k, reads[k]
  return bestg1


#################
###  CLASSES  ###
#################

class A_Bruijn_Graph():
  def __init__(self, c_reads):
    # Input:
    #   Reads stored as kt-mers and distances b/w kt-mers
    # Output:
    #   Builds a graph of nodes. Inits head node and allnodes
    # This graph assembly glues nodes together simultaneously.
    self.headNode = Node('head')
    self.allnodes = dict()      # Keys = kmers, Values = Nodes

    for j in range(len(c_reads)):
      read = c_reads[j]
      # print j, read
      ktmer = read[0]
      dist = read[1]
      current = self.headNode
      for i in range(len(ktmer)):
        if ktmer[i] in self.allnodes.keys():
          next = self.allnodes[ktmer[i]]
        else:
          next = Node(ktmer[i])
          self.allnodes[ktmer[i]] = next
        current.addOutEdge(dist[i], next)
        next.addInEdge(dist[i], current)
        next.reads.append(j)
        next.pos.append(sum(read[1][:i+1]))
        # print read[1][:i+1]
        # print next.ktmer, next.reads, next.pos
        current = next

  def split_large_nodes(self, reads_file):
    large_threshold = 5
    progress = 0
    for n in self.allnodes.keys():
      progress += 1
      if progress % 10000 == 0:
        print progress, datetime.datetime.now()
      node = self.allnodes[n]
      if len(node.reads) > large_threshold:
        reads = node.get_reads(reads_file)
        
        # Convert fasta to h, r list format
        h = []
        r = []
        for line in reads.splitlines():
          if line[0] == '>':
            if line not in h:
              h.append(line)
          else:
            if line not in r:
              r.append(line)
        h_master = copy.copy(h)

        # Find clusters
        clusters = []
        print len(r)
        while len(r) > 0:
          r1 = r[0]
          _k = 15
          cutoff = 50
          if len(r) > 1:
            headers = kmer_matching(r1, h[1:], r[1:], _k, cutoff)
          else:
            headers = []
          headers.append(h[0])
          clusters.append(headers)
          for item in headers:
            del r[h.index(item)]
            del h[h.index(item)]
          print len(r)

        print 'len clusters', len(clusters)
        print clusters

        if len(clusters) == 1:
          continue

        # Make new nodes for all but one cluster
        for i in range(len(clusters) - 1):
          print 'making new nodes:', i
          read_indices = [node.reads[h_master.index(s)] for s in clusters[i]]
          node_name = node.ktmer + '_' + str(i + 1)
          new_node = Node(node_name)
          for neighbor in copy.copy(node.outnodes):
            if bool(set(neighbor.reads) & set(read_indices)):
              e = node.get_out_edge(neighbor)
              new_node.addOutEdge(e, neighbor)
              neighbor.addInEdge(e, new_node)
              neighbor.removeInEdge(e, node)
              node.removeOutEdge(e, neighbor)
          for neighbor in copy.copy(node.innodes):
            if bool(set(neighbor.reads) & set(read_indices)):
              e = node.get_in_edge(neighbor)
              new_node.addInEdge(e, neighbor)
              neighbor.addOutEdge(e, new_node)
              neighbor.removeOutEdge(e, node)
              node.removeInEdge(e, neighbor)
          for item in read_indices:
            pos = node.remove_read(item)
            new_node.add_read(item, pos)
            # if neighbor.reads shares an element with current cluster,
            # remove neighbor from node.outnodes,
            # remove node from neighbor.innodes,
            # add neighbor to new_node.outnodes
            # add new_node to neighbor.innodes
            # ALL OF THE ABOVE also for edge distance
            # then repeat for all inward neighbors of node
            # also, copy over appropriate position information to new_node
            # then, delete position and read info from node
            # -- check if i'm missing anything else -- 
          self.allnodes[new_node.ktmer] = new_node
    return


  def joinInOutDeg1Edges(self, _k):
    # Joins two nodes linked by exactly 1 edge of length 1
    # Buggy, recommend not to use

    numJoined = 0
    traversed = set()
    nodes = set()
    for n in self.headNode.outnodes:
      nodes.add(n)
    
    while len(nodes) > 0:
      # print len(nodes)
      n = nodes.pop()
      traversed.add(n)
      # print n.ktmer, n.outDegree, n.inDegree

      # Only join nodes with indegree/outdegree of 1 and low edge distance
      if n.outDegree == 1 and n.outnodes[0].inDegree == 1 and n.outedges[0] < len(n.ktmer) / 2:
        nextNode = n.outnodes[0]
        dist = n.outedges[0]
        # print dist
        # print n.ktmer[dist:], nextNode.ktmer[:len(n.ktmer) - dist]
        # Ensure that kt-mers align. Fix edges b/w nodes
        if n.ktmer[dist:] == nextNode.ktmer[:len(n.ktmer) - dist]:
          if n.ktmer in self.allnodes:
            del self.allnodes[n.ktmer]
          if nextNode.ktmer in self.allnodes:
            del self.allnodes[nextNode.ktmer]
          # print n.ktmer, nextNode.ktmer,
          n.ktmer += nextNode.ktmer[len(n.ktmer) - dist:]
          # print n.ktmer
          numJoined += 1
          for i in range(len(nextNode.outnodes)):
            # if nextNode.outnodes[i] not in traversed:
            #   nodes.add(nextNode.outnodes[i])
            n.addOutEdge(nextNode.outedges[i] + dist, nextNode.outnodes[i])
            nextNode.outnodes[i].addInEdge(nextNode.outedges[i] + dist, n)
            nextNode.outnodes[i].removeInEdge(nextNode.outedges[i], nextNode)
          n.removeOutEdge(n.outedges[0], nextNode)
          nextNode.removeInEdge(nextNode.inedges[0], n)
          nodes.add(n)
          self.allnodes[n.ktmer] = n
        else:
          print 'ERROR:', n.ktmer, nextNode.ktmer
      else:
        for i in range(len(n.outedges)):
          if n.outnodes[i] not in traversed:
            nodes.add(n.outnodes[i])
    return numJoined

  def __str__(self):
    rep = ''
    for key in self.allnodes:
      n = self.allnodes[key]
      rep += n.ktmer + ' ' + ' '.join([str(s) for s in n.reads]) + ' ' + ' '.join([s.ktmer for s in n.outnodes]) + ' ' + ' '.join([str(s) for s in n.outedges]) + '\n'
      print len(n.outnodes), len(n.innodes)
    return rep

class Node():
  def __init__(self, ktmer):
    self.ktmer = ktmer
    self.outReadCount = list()
    self.outedges = list()
    self.outnodes = list()
    self.inedges = list()
    self.innodes = list()
    self.outDegree = 0
    self.inDegree = 0
    self.reads = list()
    self.pos = list()

  def add_read(self, read, position):
    self.reads.append(read)
    self.pos.append(position)

  def remove_read(self, target):
    # reads are stored as numbers
    if target in self.reads:
      pos_in_read = self.pos[self.reads.index(target)]
      del self.pos[self.reads.index(target)]
      del self.reads[self.reads.index(target)]
      return pos_in_read
    else:
      print 'tried to find', target, 'in', self.reads
      return None

  def get_reads(self, reads_file):
    fasta = ''
    get_line = False    
    curr_index = -1
    with open(reads_file) as f:
      for i, line in enumerate(f):
        if get_line:
          get_line = False
          fasta += line
        if line[0] == '>':
          curr_index += 1
        if curr_index in self.reads:
          fasta += line
          get_line = True
    return fasta

  def get_in_edge(self, node):
    return self.inedges[self.innodes.index(node)]

  def get_out_edge(self, node):
    return self.outedges[self.outnodes.index(node)]

  def addOutEdge(self, edge, node):
    for i in range(len(self.outedges)):
      if node == self.outnodes[i]:
        if edge == self.outedges[i]:
          self.outReadCount[i] += 1
          return
    self.outedges.append(edge)
    self.outnodes.append(node)
    self.outReadCount.append(1)
    self.outDegree += 1

  def removeOutEdge(self, edge, node):

    if node in self.outnodes:
      for i in range(len(self.outnodes)):
        if self.outnodes[i] == node and self.outedges[i] == edge:
          del self.outnodes[i]
          del self.outedges[i]
          del self.outReadCount[i]
          self.outDegree -= 1
          return True
    return False

  def getAllOutEdges(self):
    return (self.outnodes, self.outedges)

  def addInEdge(self, edge, node):
    # Do not increase indegree for head
    if node.ktmer == 'head':
      return

    for i in range(len(self.inedges)):
      if node == self.innodes[i]:
        if edge == self.inedges[i]:
          return    
    self.inedges.append(edge)
    self.innodes.append(node)
    self.inDegree += 1

  def removeInEdge(self, edge, node):
    if node in self.innodes:
      for i in range(len(self.innodes)):
        if self.innodes[i] == node and self.inedges[i] == edge:
          del self.innodes[i]
          del self.inedges[i]
          self.inDegree -= 1
          return

  def getAllInEdges(self, node):
    return (self.innodes, self.inedges)

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start