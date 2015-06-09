# Generate regular 1-deg nhoods trimmed for Yu (6/8/15)

import sys, string, datetime, random, copy, os, commands, fnmatch
import numpy as np
import itec4
from collections import defaultdict

def main():
  prior = '/home/yu/max/research/'
  reads_file = prior + 'data/reads.20k.rc.fasta'
  creads_file = prior + 'data/temp_creads.outrx_27_6_rc_v2.out'
  ktmer_headers_file = prior + 'data/temp_ktmer_headersrx_27_6_rc_v2.out'

  hr, rr = rf.read_fasta(reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  creads = build_creads_dict(creads_file, hr, rr)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]

  out_fold = '/home/yu/max/research/6.8.15_nhoods/'
  if not os.path.exists(out_fold):
    os.makedirs(out_fold)
    
  for i in [s for s in range(len(hr)) if s % 2 == 1]:
    print i
    header = hr[i]
    nhood = get_special_1_deg_nhood(header, creads, headers, hr)

    base_file = str(i) + '_base.fasta'
    hood_file = str(i) + '_hood.fasta'
    with open(base_file, 'w') as f:
      f.write(header + '\n' + rr[i])
    with open(hood_file, 'w') as f:
        f.write('\n'.join(nhood))

    commands.getstatusoutput('mv ' + base_file + ' ' + out_fold)
    commands.getstatusoutput('mv ' + hood_file + ' ' + out_fold)
  return

def get_special_1_deg_nhood(header, creads, headers, hr, rr, n_range = []):
  # Gets the special 1-deg nhood

  def convert_pos_range_to_indices(n_range, cread):
    total = 0
    indices = [-1, -1]
    for i in range(len(cread)):
      if i % 2 == 0:
        total += int(cread[i])
        if total == n_range[0] and indices[0] == -1:
          indices[0] = i
        if total == n_range[1] and indices[1] == -1:
          indices[1] = i
          break
    return indices

  def len_read(cread):
    return sum([int(cread[s]) for s in range(len(cread)) if s % 2 == 0])

  def get_pos_in_read(kmer, cread):
    if kmer in cread:
      return sum([int(cread[s]) for s in range(cread.index(kmer)) if s % 2 == 0])
    else:
      print 'error:', kmer, 'not in', cread

  # NEW FASTER? code for get_1_deg_nhood(...), test 5/25/15
  positions = defaultdict(list)    # Key = header, Val = list of ktmers and their total positions
  if header not in creads or len(creads[header]) == 1:
    print header
    return '', None
  curr_dist = 0
  master_range = []
  if len(n_range) == 0:
    master_range = range(len(creads[header]))
  else:
    indices = convert_pos_range_to_indices(n_range, creads[header])
    master_range = range(indices[0], indices[1])
  for i in master_range:
    if i % 2 == 1:
      ktmer = creads[header][i]
      for h in headers[ktmer]:
        positions[h].append(ktmer)
        positions[h].append(str(curr_dist))
    else:
      curr_dist += int(creads[header][i])

  nhood = []
  for k in positions.keys():
    curr_list = positions[k]
    dists = []
    for i in range(len(curr_list)):
      if i % 2 == 1:
        if len(dists) == 0:
          dists.append(int(curr_list[i]))
        else:
          dists.append(int(curr_list[i]) - int(curr_list[i - 2]))
    dists = dists[1:]
    # print dists

    # Stats for Yu (6/8/15)
    first_ktmer = curr_list[0]
    last_ktmer = curr_list[-2]
    first_pos_b = get_pos_in_read(first_ktmer, header)
    last_pos_b = get_pos_in_read(last_ktmer, header)
    first_pos_n = get_pos_in_read(first_ktmer, k)
    last_pos_n = get_pos_in_read(last_ktmer, k)
    dist_last_b = len_read(creads[header]) - last_pos_b
    dist_last_n = len_read(creads[k]) - last_pos_n

    if first_pos_n > first_pos_b:
      start_pos = first_pos_n - first_pos_b
      first_term = first_pos_b
    else:
      start_pos = 0
      first_term = first_pos_n

    if dist_last_n > dist_last_b:
      end_pos = len_read(creads[k]) - (dist_last_n - dist_last_b)
      last_term = dist_last_b
    else:
      end_pos = len_read(creads[k])
      last_term = dist_last_n

    nhood += [k + '_' + '_'.join([str(s) for s in [first_term, first_pos_b, last_pos_b, last_term]]), rr[hr.index(k)][start_pos : end_pos]]

  return nhood

    
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start