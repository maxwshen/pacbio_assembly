# Treats reads as the basis for error correction and extension, rather than kt-mers as in iterative_ec3.py
#
# Related files:
# combine_ec_contigs.py : Combines the read sets produced by this code into contigs

import sys, string, datetime, random, copy, os, commands, fnmatch
import numpy as np
from collections import defaultdict

import read_fasta as rf
import find_read
import kmer_matching

global temp_sig
temp_sig = str(datetime.datetime.now()).split()[1]
contigs_fold = '/home/mshen/research/contigs16/'  
overlap_accuracy_cutoff = 75    # .
overlap_length_cutoff = 300     # .
num_attempts = 2                # Number of times to try nhood extension.
support_cutoff = 70             # CANDIDATE: Required pct accuracy for support to count
support_ratio = 0.6             # CANDIDATE: Required support for a chosen read from other candidates
limit_km_times_total = 4        # How many times to attempt k-mer matching extension per direction
km_k = 15                       # .
km_cutoff = 10                  # .
support_dist_cutoff = 100000    # CONSENSUS: Bp. length, acceptable support distance from end of consensus
support_t = 3                   # CONSENSUS: Req. # reads to support a position to determine farthest support
blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
blasr_options = '-bestn 1 -m 1'   # Concise output
e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
ec_prefix = '3X_'

def main():
  reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  creads_file = '/home/mshen/research/data/22.4_creads.out'
  ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  ec_tool = '/home/mshen/research/bin/error_correction_3X_0112.sh'
  parallel_prefix = sys.argv[1]
  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', ktmer_headers_file, '\nEC Tool:', ec_tool


  # Actions
  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix)
  # ktmer_reads_pct_overlap(ktmer_headers_file, reads_file)
  # combine_contigs(contigs_fold)
  # check_contigs(contigs_fold, reads_file)

def check_contigs(contigs_fold, reads_file):
  # Aligns component reads to combined contigs in an effort to find jumps.
  # 1/5/15: Doesn't work very well.
  hr, rr = rf.read_fasta(reads_file)
  res_file = contigs_fold + 'contig_0results.fasta'
  comb_file = contigs_fold + 'contig_0_combined.fasta'
  temp_file = 'temp_cc_' + temp_sig + '.fasta'
  with open(res_file) as f:
    reads = [s.split()[0] for s in f.readlines()]
  for r in reads:
    mod_r = '/'.join(r.split('/')[:-1]).replace('START', '')
    print mod_r
    with open(temp_file, 'w') as f:
      f.write(mod_r + '\n' + rr[hr.index(mod_r)])
    status = commands.getstatusoutput(blasr_exe + ' ' + temp_file + ' ' + comb_file + ' ' + blasr_options)[1]
    print status


def combine_contigs(contigs_fold):
  acc_cutoff = 0
  dist_to_end = 100
  temp_base = 'temp_contigbase_' + temp_sig + '.fasta'
  temp_try = 'temp_contigtry_' + temp_sig + '.fasta'
  for fn in os.listdir(contigs_fold):
    if fnmatch.fnmatch(fn, '*[0-9].fasta'):
      print fn                              # TESTING
      hs, rs = rf.read_fasta(contigs_fold + fn)
      bases = [rs[0]]
      for i in range(1, len(rs)):
        new_bases = []
        num_fails = 0
        print len(bases)
        for base in bases:
          new_base = base
          with open(temp_base, 'w') as f:
            f.write('>base\n' + base)
          with open(temp_try, 'w') as f:
            f.write('>try\n' + rs[i])
          status = commands.getstatusoutput(blasr_exe + ' ' + temp_try + ' ' + temp_base + ' ' + blasr_options)[1]
          if len(status) == 0:
            print 'FAILED BLASR ALIGNMENT'
            num_fails += 1
            new_bases.append(new_base)
            continue
          else:
            print status                      # TESTING
            acc = float(status.split()[5])
            beg_align_r1 = int(status.split()[6])
            end_align_r1 = int(status.split()[7])
            total_len_r1 = int(status.split()[8])
            end_pos_r1 = total_len_r1 - end_align_r1
            beg_align_r2 = int(status.split()[9])
            end_align_r2 = int(status.split()[10])
            total_len_r2 = int(status.split()[11])
            end_pos_r2 = total_len_r2 - end_align_r2

            if end_pos_r1 < dist_to_end: # and beg_align_r2 < dist_to_end:
              if acc > acc_cutoff:
                new_base = base[: end_align_r1]
                new_base += rs[i][end_align_r2 : ]
            new_bases.append(new_base)
        if num_fails == len(bases):
          new_bases.append(rs[i])
        bases = new_bases

      best_base = sorted(bases, key = len, reverse = True)[0]
      out_file = fn.split('.')[0] + '_combined.fasta'
      with open(contigs_fold + out_file, 'w') as f:
        f.write('>' + fn + '\n' + best_base)
      status = commands.getstatusoutput(blasr_exe + ' ' + contigs_fold + out_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]
      print 'BEST:'
      print status


def ktmer_reads_pct_overlap(ktmer_headers_file, reads_file):
  # Finds read clusters (aligned to the genome) for all kt-mers
  def within(beg1, end1, beg2, end2):
    # Test if 2 is in 1
    if beg1 < beg2 < end1:
      return True
    if beg1 < end2 < end1:
      return True
    return False

  def expand(beg1, end1, beg2, end2):
    newbeg1 = beg1
    newend1 = end1
    if beg2 < beg1:
      newbeg1 = beg2
    if end2 > end1:
      newend1 = end2
    return [newbeg1, newend1]

  def within_all(clusters):
    for i in range(len(clusters)):
      for j in range(len(clusters)):
        if i != j:
          if within(clusters[i][0], clusters[i][1], clusters[j][0], clusters[j][1]):
            return i, j
    return None

  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)

  for kt in headers.keys():
    print kt
    clusters = []
    for h in headers[kt]:
      tempfile = 'temp' + temp_sig + '.fasta'
      with open(tempfile, 'w') as f:
        f.write('>1\n' + rr[hr.index(h)]) 
      status = commands.getstatusoutput(blasr_exe + ' ' + tempfile + ' ' + e_coli_genome + ' ' + blasr_options)[1]
      if len(status) > 0:
        beg = int(status.split()[6])
        end = int(status.split()[7])
        found = False
        for c in clusters:
          if within(c[0], c[1], beg, end):
            found = True
            new = expand(c[0], c[1], beg, end)
            c[0] = new[0]
            c[1] = new[1]
            c[2] += 1
        if not found:
          clusters.append([beg, end, 1])
        # print clusters, beg, end

    # Do final clustering of clusters
    while True:
      found = False
      if within_all(clusters) is not None:
        found = True
        i, j = within_all(clusters)
        new = expand(clusters[i][0], clusters[i][1], clusters[j][0], clusters[j][1])
        clusters[i][0] = new[0]
        clusters[i][1] = new[1]
        clusters[i][2] += clusters[j][2]
        clusters.remove(clusters[j])
      if not found:
        break

    # Currently just prints genomic alignment of clusters.
    for c in clusters:
      print ' '.join([str(s) for s in c])

    # May want to write a sister function that doesn't rely on the genome


def verify_against_candidates(h, candidate_headers, hr, rr):
  # I think this can be improved to detect jumps better 
  dist_to_end = 100
  if len(candidate_headers) == 0:
    return True

  temp_file = 'temp_h_' + temp_sig + '.fasta'
  with open(temp_file, 'w') as f:
    f.write(h + '\n' + rr[hr.index(h)])

  support = 0
  temp_chfile = 'temp_ch_' + temp_sig + '.fasta'
  for ch in candidate_headers:
    with open(temp_chfile, 'w') as f:
      f.write(ch + '\n' + rr[hr.index(ch)])
    status = commands.getstatusoutput(blasr_exe + ' ' + temp_file +' ' + temp_chfile + ' ' + blasr_options)[1]
    print status                        # TESTING
    if len(status) != 0:
      acc = float(status.split()[5])
      beg_align_r1 = int(status.split()[6])
      end_align_r1 = int(status.split()[7])
      total_len_r1 = int(status.split()[8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[9])
      end_align_r2 = int(status.split()[10])
      total_len_r2 = int(status.split()[11])
      end_pos_r2 = total_len_r2 - end_align_r2
      length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2     # Avg alignment length

      is_support = False
      if acc > support_cutoff and length > support_len:
        is_support = True
      if end_pos_r1 < dist_to_end and beg_align_r2 < dist_to_end:
        is_support = True
      if end_pos_r2 < dist_to_end and beg_align_r1 < dist_to_end:
        is_support = True
      if end_pos_r1 < dist_to_end and beg_align_r1 < dist_to_end:
        is_support = True
      if end_pos_r2 < dist_to_end and beg_align_r2 < dist_to_end:
        is_support = True

      if is_support:
          support += 1

  support_pct = float(support) / float(len(candidate_headers))
  print 'support found:', support_pct         # TESTING
  return support_pct >= support_ratio


def verify_against_candidates_longest(candidate_headers, hr, rr):
  # Improved version
  # Take the longest read, find support against that - no long tails
  dist_to_end = 300
  if len(candidate_headers) == 0:
    return True

  cand_r = [rr[hr.index(s)] for s in candidate_headers]
  cand_r.sort(key = len, reverse = True)
  base = cand_r[0]

  temp_base = 'temp_bas_' + temp_sig + '.fasta'
  with open(temp_base, 'w') as f:
    f.write('>' + hr[rr.index(base)] + '\n' + base)

  print commands.getstatusoutput(blasr_exe + ' ' + temp_base +' ' + e_coli_genome + ' ' + blasr_options)[1]

  support = 0
  temp_chfile = 'temp_ch_' + temp_sig + '.fasta'
  for ch in cand_r[1:]:
    with open(temp_chfile, 'w') as f:
      f.write('>' + hr[rr.index(ch)] + '\n' + ch)
    status = commands.getstatusoutput(blasr_exe + ' ' + temp_base +' ' + temp_chfile + ' ' + blasr_options)[1]
    print status                        # TESTING
    print commands.getstatusoutput(blasr_exe + ' ' + temp_chfile +' ' + e_coli_genome + ' ' + blasr_options)[1]
    if len(status) != 0:
      acc = float(status.split()[5])
      beg_align_r1 = int(status.split()[6])
      end_align_r1 = int(status.split()[7])
      total_len_r1 = int(status.split()[8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[9])
      end_align_r2 = int(status.split()[10])
      total_len_r2 = int(status.split()[11])
      end_pos_r2 = total_len_r2 - end_align_r2
      length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2     # Avg alignment length

      is_support = False
      if acc > support_cutoff:
        if beg_align_r1 > dist_to_end and beg_align_r2 < dist_to_end or beg_align_r1 < dist_to_end and beg_align_r2 > dist_to_end or beg_align_r1 < dist_to_end and beg_align_r2 < dist_to_end:
          if end_pos_r1 > dist_to_end and end_pos_r2 < dist_to_end or end_pos_r1 < dist_to_end and end_pos_r2 > dist_to_end or end_pos_r1 < dist_to_end and end_pos_r2 < dist_to_end: 
            is_support = True
            support += 1
            print '+1'

  support_pct = float(support) / float(len(candidate_headers) - 1)
  print 'support found:', support_pct         # TESTING
  return support_pct >= support_ratio


def filter_ktmers(ktmers, creads, headers):
  # Filters out kt-mers that are in the same 1-deg nhood
  # Produces about 3000 from 240k 22,4-mers, evenly spaced
  # Max dist 12k, average 1.5k
  new_ktmers = []
  excluded = set()
  for kt in ktmers:
    if kt in excluded:
      continue
    new_ktmers.append(kt)
    for n in find_neighboring_ktmers(kt, headers, creads):
      excluded.add(n)
  return new_ktmers

  # genome = '/home/mshen/research/data/e_coli_genome.fasta'
  # gh, gr = rf.read_fasta(genome)
  # for kt in ktmers:
  #   if kt in gr[0]:
  #     print gr[0].index(kt)


def ktmers_from_genome(ktmers, min_bp, max_bp):
  hg, rg = rf.read_fasta(e_coli_genome)
  rg = rg[0]

  new_kt = []
  _k = len(ktmers[0])
  for i in range(min_bp, max_bp - _k + 1):
    if rg[i : i + _k] in ktmers:
      new_kt.append(rg[i : i + _k])
  return new_kt


def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix):
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'

  contigs = []

  curr_min_pos = 0
  curr_max_pos = 5000000
  if parallel_prefix == '000':
    curr_max_pos = 500000
  if parallel_prefix == '050':
    curr_min_pos = 500000
    curr_max_pos = 1000000
  if parallel_prefix == '100':
    curr_min_pos = 1000000
    curr_max_pos = 1500000
  if parallel_prefix == '150':
    curr_min_pos = 1500000
    curr_max_pos = 2000000
  if parallel_prefix == '200':
    curr_min_pos = 2000000
    curr_max_pos = 2500000
  if parallel_prefix == '250':
    curr_min_pos = 2500000
    curr_max_pos = 3000000
  if parallel_prefix == '300':
    curr_min_pos = 3000000
    curr_max_pos = 3500000
  if parallel_prefix == '350':
    curr_min_pos = 3500000
    curr_max_pos = 4000000
  if parallel_prefix == '400':
    curr_min_pos = 4000000
    curr_max_pos = 4500000
  if parallel_prefix == '450':
    curr_min_pos = 4500000
    curr_max_pos = 5000000

  # min_bp = 106713
  # max_bp = 108843
  # ktmers = ktmers_from_genome(ktmers, min_bp, max_bp)   # testing
  covered_range = []        # testing, stores a list of consensus positions so we don't overlap
  ktmers = filter_ktmers(ktmers, creads, headers)
  print 'After filtering,', len(ktmers), 'kt-mers remain.'

  # num_contig_attempts = 200                   # testing
  num_contig_attempts = len(ktmers)
  for m in range(num_contig_attempts):
    print '\n' + str(datetime.datetime.now())
    curr_ktmer = ktmers[m]

    h = get_read_with_most_neighbors(curr_ktmer, headers, creads)
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/25524/0_5702/0_5702'  # jump ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/66451/0_6519/0_6519'    # normal ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/22681/0_4859/0_4859'  # 4 iterations  before jump ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/1326/0_5814/0_5814' # a little farther from jump ex
    print 'STARTING HEADER:\n', h
    curr_contig = [error_correct(ec_tool, h, headers, creads, hr, rr)]
    print 'STARTING AT',                        # testing
    pos = find_genomic_position(curr_contig[0])       # testing
    if pos > curr_max_pos or pos < curr_min_pos:                   # testing
      continue                                  # testing
    for cr in covered_range:
      if abs(pos - cr) < 1500:
        # don't start here, if it's within 1500bp of a consensus we already have 
        continue 
    curr_contig_headers = [h + 'START']
    master_h = h
    master_traversed_headers = [h]


    # MAIN LOOP
    for direction in ['right', 'left']:
      counter = 0
      limit_km_times = limit_km_times_total
      h = master_h
      while True:
        # Break condition: Current header doesn't change, meaning we couldn't find any extension candidates
        counter += 1
        print 'iteration', counter, direction
        old_h = h
        temp_traversed_headers = []
        if limit_km_times > 0:
          num_attempts_temp = num_attempts + 1
        else:
          num_attempts_temp = num_attempts
        for i in range(num_attempts_temp):
          print 'Attempt', i                                            # testing
          km = False
          km_early_out = False
          traversed_headers = master_traversed_headers + temp_traversed_headers

          # Grab candidates via nhood extension or kmer matching
          if i < num_attempts:
            # Try nhood extension num_attempt times, then try kmer matching if still no extension
            possible_heads = extend_n(h, headers, creads, traversed_headers, direction, hr, rr)
            print len(possible_heads), 'candidates for extension'
            if len(possible_heads) == 0:
              print 'could not extend further'
              continue
          else:
            if km_early_out:
              possible_heads = []
              break
            print 'Trying k-mer matching'
            limit_km_times -= 1
            km = True
            if direction == 'right':
              possible_heads = kmer_matching.kmer_matching(curr_contig[-1], reads_file, km_k, km_cutoff, file_bool = False)
            if direction == 'left':
              possible_heads = kmer_matching.kmer_matching(curr_contig[0], reads_file, km_k, km_cutoff, file_bool = False)
            if h in possible_heads:
              possible_heads.remove(h)

          # Filter candidates by overlapping test
          # If we find farthest_support is far, redo overlapping with trimmed consensus
          criteria = dict()   # Key = header, Val = some order-able criteria (ex: length)
          for i2 in range(2):
            good_candidates = []
            farthest_support = []
            for head in possible_heads:
              candidate_read = rr[hr.index(head)]
              # if len(possible_heads) < 50:                            # testing
                # print head,                                           # testing
                # find_genomic_position(candidate_read)       # testing

              # candidate_read = error_correct(ec_tool, head, headers, creads, hr, rr)    # testing
              # print 'consensus:',                         # testing
              # find_genomic_position(candidate_read)       # testing
              # print' original:',                          # testing
              # find_genomic_position(rr[hr.index(head)])   # testing

              overlaps = False
              if direction == 'right' and test_overlap(head, candidate_read, curr_contig[-1], direction, farthest_support, criteria, relaxed = False):
                # option: relaxed = km
                overlaps = True
              if direction == 'left' and test_overlap(head, curr_contig[0], candidate_read, direction, farthest_support, criteria, relaxed = False):
                overlaps = True
              if overlaps:
                good_candidates.append(head)
              # print 'Overlap:', overlaps                  # testing

            if len(farthest_support) == 0:
              break
            support_pos = min(support_t - 1, len(farthest_support) - 1)
            farthest_support.sort()
            best_support = farthest_support[support_pos]
            if best_support < support_dist_cutoff:
              break
            else:
              print 'Trimmed consensus since support dist =', farthest_support[support_pos]
              if direction == 'right':
                curr_contig[-1] = curr_contig[-1][ : -farthest_support[support_pos]]
              if direction == 'left':
                curr_contig[0] = curr_contig[0][farthest_support[support_pos] : ]
          print len(good_candidates), 'reads passed overlap filter'

          # for gc in good_candidates:                                # testing
            # find_genomic_position(rr[hr.index(gc)])       # testing

          # Filter candidates by their support for each other
          filtered_good_candidates = []
          if not km:
            filtered_good_candidates = good_candidates
            # if verify_against_candidates_longest(good_candidates, hr, rr):
              # filtered_good_candidates = good_candidates
            # for gc in good_candidates:
              # temp_candidates = copy.copy(good_candidates)
              # temp_candidates.remove(gc)
              # if verify_against_candidates(gc, temp_candidates, hr, rr):
                # filtered_good_candidates.append(gc)
          if km:
            filtered_good_candidates = good_candidates
          print 'Filtered', len(good_candidates) - len(filtered_good_candidates), 'reads by support filter'
          if len(filtered_good_candidates) == 0:
            print 'No reads passed filtering step against other candidates'
            for p in possible_heads:
              temp_traversed_headers.append(p)
            if km:
              km_early_out = True
            continue

          # Sort candidates by some criteria
          filtered_good_candidates.sort(key = lambda d: criteria[d], reverse = True)
          random.shuffle(filtered_good_candidates)    # testing

          # Once we choose a particular candidate
          consensus_temp = ''
          for i in range(len(filtered_good_candidates)):
            best_head = filtered_good_candidates[i]
            h = best_head
            consensus_temp = error_correct(ec_tool, h, headers, creads, hr, rr)
            if len(consensus_temp) != 0 and consensus_temp not in curr_contig:
              break  
          if len(consensus_temp) == 0:
            print 'COULD NOT ERROR CORRECT ANY FILTERED GOOD CANDIDATES'
            h = old_h
          else:
            print 'New header:', h, criteria[h]       # testing
            print 'SUCCESS!',                         # testing 
            con_pos = find_genomic_position(consensus_temp)     # testing
            covered_range.append(con_pos)
            if len(consensus_temp) == 0:
              print 'failed to error correct'
              continue
            if direction == 'right':
              curr_contig.append(consensus_temp)
              curr_contig_headers.append(h)
            if direction == 'left':
              curr_contig.insert(0, consensus_temp)
              curr_contig_headers.insert(0, h)
            master_traversed_headers.append(h)

          if h == old_h:
            # This part ensures nhood extension by preventing us from 
            # getting stuck in returning lower degree headers. This should be made more temporary.
            for p in possible_heads:
              temp_traversed_headers.append(p)
          else:
            break
        if h == old_h:
          print 'No new reads overlapped with current contig'
          break

    # ASSESS RESULTS
    print old_h
    contig = ''
    contig_file = contigs_fold + 'contig_' + parallel_prefix + str(m) + '.fasta'
    contig_result = contigs_fold + 'contig_' + parallel_prefix + str(m) + 'results.fasta'
    for j in range(len(curr_contig)):
      if curr_contig[j] != '':
        contig += '>' + curr_contig_headers[j] + '\n' + curr_contig[j] + '\n'
    with open(contig_file, 'w') as f:
      f.write(contig)
    status = commands.getstatusoutput(blasr_exe + ' ' + contig_file +' ' + e_coli_genome + ' ' + blasr_options + ' > ' + contig_result)[1]

    # Filter kt-mers
    curr_ktmer_len = len(ktmers)
    for i in range(len(curr_contig)):
      new_ktmers = []
      for kt in ktmers:
        if kt not in curr_contig[i]:
          new_ktmers.append(kt)
      ktmers = new_ktmers
    print curr_ktmer_len - len(ktmers), ' kt-mers filtered from consensus'


def find_genomic_position(read):
  temp_file = 'temp_read' + temp_sig + '.fasta'
  with open(temp_file, 'w') as f:
    f.write('>1\n' + read)
  status = commands.getstatusoutput(blasr_exe + ' ' + temp_file +' ' + e_coli_genome + ' ' + blasr_options)[1]

  if len(status) > 0:
    beg = int(status.split()[6])
    end = int(status.split()[7])
    print '\taligned to:', beg, end
    return (beg + end) / 2
  else:
    print '\tFAILED ALIGNMENT'
    return -1
  


def test_overlap(head1, seq1, seq2, direction, farthest_support, criteria, relaxed = False):
  # Tests that seq1 is after seq2
  # farthest_support is a list that will contains distances 
  # from the end (depending on direction) of the current read
  dist_from_end = 100
  acc_cutoff = overlap_accuracy_cutoff
  len_cutoff = overlap_length_cutoff

  temps1 = 'temp_seq1' + temp_sig + '.fasta'
  temps2 = 'temp_seq2' + temp_sig + '.fasta'
  with open(temps1, 'w') as f:
    f.write('>1\n' + seq1)
  with open(temps2, 'w') as f:
    f.write('>2\n' + seq2)

  status = commands.getstatusoutput(blasr_exe + ' ' + temps1 + ' ' + temps2 + ' ' + blasr_options)[1]
  if len(status.strip()) == 0:
    return False
  # print status                    # TESTING
  accuracy = float(status.split()[5])
  beg_align_r1 = int(status.split()[6])
  end_align_r1 = int(status.split()[7])
  total_len_r1 = int(status.split()[8])
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[9])
  end_align_r2 = int(status.split()[10])
  total_len_r2 = int(status.split()[11])
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2          # Average alignment length

  # update criteria
  criteria[head1] = length

  # Update farthest support, the distance to the end of the consensus that has support from 1deg nhood 
  if direction == 'right':
    farthest_support.append(end_pos_r1)
  if direction == 'left':
    farthest_support.append(beg_align_r2)

  if not relaxed:
    if direction == 'right':
      # if accuracy >= acc_cutoff and length > len_cutoff and end_pos_r1 < dist_from_end and end_pos_r2 > end_pos_r1 and beg_align_r2 < dist_from_end:
        # print status                    # TESTING
      return accuracy >= acc_cutoff and length > len_cutoff and end_pos_r1 < dist_from_end and end_pos_r2 > end_pos_r1 and beg_align_r2 < dist_from_end
    if direction == 'left':
      return accuracy >= acc_cutoff and length > len_cutoff and beg_align_r2 < dist_from_end and beg_align_r1 > beg_align_r2 and end_pos_r1 < dist_from_end
  else:
    return accuracy >= acc_cutoff and length > len_cutoff


def extend_n(header, headers, creads, traversed_headers, direction, hr, rr):
  leniency = 100    # If 1-degree nhood fails, we accept a read that extends beyond border - leniency
  backtrack_limit = 1000       # Don't backtrack too far
  num_kmers_cutoff = 500       # If we start considering this many kt-mers in nhood extension, stop.
  
  dist_to_end = get_dist_to_end(header, creads, direction)
  ktmers = dist_to_end.keys()
  print len(ktmers), 'ktmers in 1-deg nhood'

  reads = []
  num_neighbors = []
  for k in ktmers:
    next_read = find_extending_read(k, headers, hr, rr)
    if len(next_read) != 0:
      accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
      reads += accepted
      for nr in next_read:
        num_neighbors.append(len(creads[nr]) / 2 - 1)
  if len(reads) != 0:
    return [x for (y, x) in sorted(zip(num_neighbors, reads), reverse = True)]

  # Try finding a read that comes close to passing the current read
  for k in dist_to_end.keys():
    dist_to_end[k] -= leniency

  print 'trying n-degree nhood extension'
  traversed_ktmers = set(ktmers)
  while True:
    print '\tkmers considered:', len(traversed_ktmers)
    new_ktmers = defaultdict(list)   # Key = new-ktmer, Val = [old ktmer that is connected]
    if len(ktmers) > num_kmers_cutoff:
      print 'Stopping - too many kt-mers'
      return ''
    for kt in ktmers:
      for kn in find_neighboring_ktmers(kt, headers, creads):
        if kn not in traversed_ktmers:
          if kn not in new_ktmers.keys() or kt not in new_ktmers[kn]:
            new_ktmers[kn].append(kt)
    if len(new_ktmers.keys()) == 0:
      print 'None found'
      return ''

    print 'found', len(new_ktmers.keys()), 'ktmers from', len(ktmers), 'ktmers'
    new_dist_to_end = dict()
    reads = []
    num_neighbors = []
    for kt in new_ktmers.keys():
      k = new_ktmers[kt][0]
      dist = dist_bw_ktmers(kt, k, headers, creads)
      if dist == None:
        continue
      if direction == 'right':
        new_dist_to_end[kt] = dist_to_end[k] + dist
      if direction == 'left':
        new_dist_to_end[kt] = dist_to_end[k] - dist
      next_read = find_extending_read(k, headers, hr, rr)
      if len(next_read) != 0:
        accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
        reads += accepted
        for nr in next_read:
          num_neighbors.append(len(creads[nr]) / 2 - 1)
    if len(reads) != 0:
      return [x for (y, x) in sorted(zip(num_neighbors, reads), reverse = True)]

    # Filter out those who are too far
    for k in ktmers:
      traversed_ktmers.add(k)

    # print new_dist_to_end
    ktmers = []
    for kt in new_ktmers.keys():
      if new_dist_to_end[kt] < backtrack_limit and kt not in traversed_ktmers:
        ktmers.append(kt)

    dist_to_end = new_dist_to_end


def find_neighboring_ktmers(ktmer, headers, creads):
  # Speed bottleneck 12/25/14
  # But storing dict of all neighbors in memory is slower
  neighbors = []
  for h in headers[ktmer]:
    for i in range(len(creads[h])):
      if i % 2 == 1:
        neighbors.append(creads[h][i])
  return neighbors


def get_dist_to_end(header, creads, direction):
  dist_to_end = dict()    # Key = ktmer in read, Val = distance to end of read
  if direction == 'right':
    for i in range(len(creads[header])):
      if i % 2 == 1:
        dist_to_end[creads[header][i]] = 0
      elif i > 0:
        for k in dist_to_end:
          dist_to_end[k] += int(creads[header][i])
  if direction == 'left':
    curr = 0
    for i in range(len(creads[header])):
      if i % 2 == 1:
        dist_to_end[creads[header][i]] = curr
      else:
        curr += int(creads[header][i])
  return dist_to_end


def dist_bw_ktmers(kt1, kt2, headers, creads):
  # Finds distance between kt1 and kt2 in a read. If kt2 is before kt1, dist is negative
  found = False
  for h in headers[kt1]:
    if kt2 in creads[h]:
      header = h
      found = True

  if not found:
    return None

  dist = 0
  grab = False
  for i in range(len(creads[header])):
    if creads[header][i] == kt1 or creads[header][i] == kt2:
      if grab == True:
        grab = False
        break
      if grab == False:
        grab = True
    if grab and i % 2 == 0:
      dist += int(creads[header][i])
  if creads[header].index(kt1) > creads[header].index(kt2):
    dist *= -1
  return dist


def find_extending_read(ktmer, headers, hr, rr):
  # If a read in ktmer extends past dist, return it
  # Modified 12/31/14: Do not check if read extends past dist, instead
  # check this to the last error corrected portion.
  # This is because sometimes the EC consensus can be much shorter than
  # the actual read, and we just want a read that extends past the EC consensus,
  # not the current read we perform 1-deg nhood extension on.
  valid = []
  # if verify_against_candidates(headers[ktmer][0], headers[ktmer][1:], hr, rr):
    # valid = headers[ktmer]
  valid = headers[ktmer]  
  return valid


def error_correct(ec_tool, header, headers, creads, hr, rr, temp_sig_out = None):
  if temp_sig_out is not None:
    temp_sig = temp_sig_out
  else:
    global temp_sig
  reads = []
  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    return ''
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      collected_h.add(h)

  reads = [header, rr[hr.index(header)]]
  for ch in collected_h:
    reads.append(ch)
    reads.append(rr[hr.index(ch)])
  print len(reads) / 2, 'reads used for correction'

  temp_orig_file = 'temp_orig' + temp_sig + '.fasta'
  with open(temp_orig_file, 'w') as f:
    f.write(header + '\n' + rr[hr.index(header)])

  temp_nhood_file = 'temp_nhood' + temp_sig + '.fasta'
  with open(temp_nhood_file, 'w') as f:
    f.write('\n'.join(reads))

  ec_out = ec_prefix + temp_orig_file
  status = commands.getstatusoutput(ec_tool + ' ' + temp_orig_file + ' ' + temp_nhood_file)[1]
  print status
  if 'ERROR' in status or 'No such file or directory' in status:
    print status
    return ''

  with open(ec_out, 'r') as f:  
    consensus = f.readlines()[1].strip()
  print 'consensus len:', len(consensus), 'out of', len(rr[hr.index(header)])
  return consensus

def get_read_with_most_neighbors(ktmer, headers, creads):
  best_h = ''
  best_neighbors = 0
  for h in headers[ktmer]:
    num_neighbors = len(creads[h]) / 2 - 1
    print h, num_neighbors
    if num_neighbors > best_neighbors:
      best_neighbors = num_neighbors
      best_h = h
  return best_h

def build_creads_dict(creads_file, reads_file):
  creads = defaultdict(list)   # Key = header, Val = creads 
  h, r = rf.read_fasta(reads_file)
  with open(creads_file) as f:
    for i, line in enumerate(f):
      creads[h[i]] = line.split()
  return creads

def build_headers_dict(ktmer_headers_file):
  headers = defaultdict(list)   # Key = ktmer, Val = [headers]
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      headers[words[0]] = words[1:]
  return headers

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start