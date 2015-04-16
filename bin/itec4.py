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
prior = '/home/yu/max/research/'
# contigs_fold = '/home/mshen/research/contigs_50x_1/'
# contigs_fold = '/home/max/research/contigs_50x_3/'
contigs_fold = prior + 'contigs_20kb_full_14/'
overlap_accuracy_cutoff = 75    # .
overlap_length_cutoff = 7000     # .
# overlap_length_cutoff = 300     # .
num_attempts = 1                # Number of times to try nhood extension.
support_cutoff = 70             # CANDIDATE: Required pct accuracy for support to count
support_ratio = 0.6             # CANDIDATE: Required support for a chosen read from other candidates
limit_km_times_total = 5        # How many times to attempt k-mer matching extension per direction
km_k = 15                       # .
km_cutoff = 100                 # .
support_dist_cutoff = 100000    # CONSENSUS: Bp. length, acceptable support distance from end of consensus
support_t = 3                   # CONSENSUS: Req. # reads to support a position to determine farthest support
nhood_header_limit = float('inf')         # .
nhood_it_limit = 3              # .
# blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
blasr_exe = 'blasr'
blasr_zero = 4      # 0 on debruijn, 4 on Yu's computer
blasr_zero_len = 8  # 0 on debruijn, 4 on Yu's computer
blasr_options = '-bestn 1 -m 1'   # Concise output
# e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
# e_coli_genome = '/home/max/research/data/e_coli_genome.fasta'
e_coli_genome = '/home/yu/e_coli_genome.fasta'
# ec_prefix = '3X_'
# ec_prefix = '5X_'
ec_prefix = 'C0413_'
use_ecs = False

def main():
  if not os.path.exists(contigs_fold):
    os.makedirs(contigs_fold)

  # OLD DATASET
  # reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  # creads_file = '/home/mshen/research/data/22.4_creads.out'
  # ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'

  # parallel_prefix = sys.argv[1]
  parallel_prefix = str(0)
  # cov = sys.argv[1]
  # _k = sys.argv[2]
  # _t = sys.argv[3]

  # NEW 20KB DATASET
  # reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  # reads_file = '/home/mshen/research/data/reads.20k.rc.fasta'
  # reads_file = '/home/mshen/research/data/reads.20k.' + cov + 'x.rc.fasta'
  # reads_file = '/home/max/research/data/reads.20k.' + cov + 'x.rc.fasta'
  reads_file = prior + 'data/reads.20k.rc.fasta'

  # creads_file = '/home/mshen/research/data/temp_creads.out_28_6_rc.out'
  # ktmer_headers_file = '/home/mshen/research/data/temp_ktmer_headers_28_6_rc.out'
  # creads_file = '/home/mshen/research/data/temp_creads.out' + cov + 'x_' + _k + '_' + _t + '_rc.out'
  # ktmer_headers_file = '/home/mshen/research/data/temp_ktmer_headers' + cov + 'x_' + _k + '_' + _t + '_rc.out'
  # creads_file = '/home/max/research/data/temp_creads.out' + cov + 'x_' + _k + '_' + _t + '_rc.out'
  # ktmer_headers_file = '/home/max/research/data/temp_ktmer_headers' + cov + 'x_' + _k + '_' + _t + '_rc.out'

  creads_file = prior + 'data/temp_creads.out_28_6_rc.out'
  ktmer_headers_file = prior + 'data/temp_ktmer_headers_28_6_rc.out'

  # ec_tool = '/home/lin/program/error_correction_5X_0210.sh'   
  # ec_tool = '/home/max/program/error_correction_0318.sh'      # yu's comp
  # ec_tool = '/home/yu/program/error_correction_0402.sh'
  # ec_tool = '/home/yu/program/error_correction_test.sh'
  ec_tool = '/home/yu/program/error_correction_0413.sh'

  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', \
    ktmer_headers_file, '\nEC Tool:', ec_tool, '\nContigs fold', contigs_fold


  # Actions
  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix)
  # ktmer_reads_pct_overlap(ktmer_headers_file, reads_file)
  # combine_contigs(contigs_fold)
  # check_contigs(contigs_fold, reads_file)
  # output_all_1_deg_nhoods(reads_file, creads_file, ktmer_headers_file, ec_tool, parallel_prefix)
  # contigs_results_file = '/home/mshen/research/contigs30/contig_70results.fasta'
  # output_some_1_deg_nhoods(contigs_results_file, reads_file, creads_file, ktmer_headers_file, \
    # ec_tool)
  # find_jumps_in_contigs(contigs_fold, parallel_prefix)


def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix):
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  if use_ecs:
    ecs = read_ec_from_file()
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

  # min_bp = 2119287 
  # max_bp = 2119820
  # print 'Filtering kt-mers between', min_bp, max_bp
  # ktmers = ktmers_from_genome(ktmers, min_bp, max_bp)   # testing
  covered_range = []        # testing, stores a list of consensus positions so we don't overlap
  # ktmers = filter_ktmers(ktmers, creads, headers)
  print 'After filtering,', len(ktmers), 'kt-mers remain.'

  num_contig_attempts = 10                   # testing
  # num_contig_attempts = len(ktmers)
  for m in range(num_contig_attempts):
    print '\n' + str(datetime.datetime.now())
    curr_ktmer = ktmers[m]

    h = get_read_with_most_neighbors(curr_ktmer, headers, creads)
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/4937/4611_8608/0_3997'  # farthest
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/37843/0_2510/0_2510'    # base
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/74540/0_1570/0_1570'  # most ktmers
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/25524/0_5702/0_5702'  # jump ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/66451/0_6519/0_6519'    # normal ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/22681/0_4859/0_4859'  # 4 iterations  before jump ex
    # h = '>m120114_011938_42177_c100247042550000001523002504251220_s1_p0/1326/0_5814/0_5814' # a little farther from jump ex
    print 'STARTING HEADER:\n', h
    curr_contig = [error_correct(ec_tool, h, headers, creads, hr, rr)]
    print 'STARTING AT',                        # testing
    if len(curr_contig[0]) == 0:
      continue
    pos = find_genomic_position(curr_contig[0], hr, rr, align_consensus = True)       # testing
    if pos > curr_max_pos or pos < curr_min_pos:                   # testing
      continue                                  # testing
    for cr in covered_range:
      if abs(pos - cr) < 50:
        # don't start here, if it's within 1500bp of a consensus we already have 
        continue 
    curr_contig_headers = [h + 'START']
    master_h = h
    master_traversed_headers = [h]

    curr_time = datetime.datetime.now()


    # MAIN LOOP
    # for direction in ['right', 'left']:
    for direction in ['right']:
      counter = 0
      limit_km_times = limit_km_times_total
      h = master_h
      while True:
        # Break condition: Current header doesn't change, meaning we couldn't find any extension candidates
        counter += 1
        print datetime.datetime.now() - curr_time
        print 'iteration', counter, direction
        curr_time = datetime.datetime.now()
        if counter > 500:
          break
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
            new_possible_heads = []
            for ph in possible_heads:
              new_possible_heads.append(ph.split()[0])
            possible_heads = new_possible_heads
            if h in possible_heads:
              possible_heads.remove(h)
            print 'from k-mer matching, found', len(possible_heads), 'possible reads'


          # Filter candidates by overlapping test
          # If we find farthest_support is far, redo overlapping with trimmed consensus
          criteria = dict()   # Key = header, Val = some order-able criteria (ex: length)
          for i2 in range(2):
            good_candidates = []
            farthest_support = []
            for head in possible_heads:
              if use_ecs and head in ecs:
                candidate_read = ecs[head]
              else:
                candidate_read = rr[hr.index(head)]

              # if len(possible_heads) < 50:                            # testing
                # print head,                                           # testing
                # find_genomic_position(candidate_read, hr, rr)       # testing
                # num_neighbors = len(creads[head]) / 2 - 1   # testing
                # print head, num_neighbors                   # testing

              # candidate_read = error_correct(ec_tool, head, headers, creads, hr, rr)    # testing
              # print 'consensus:',                         # testing
              # find_genomic_position(candidate_read, hr, rr, align_consensus = True))       # testing
              # print' original:',                          # testing
              # find_genomic_position(rr[hr.index(head)], hr, rr)   # testing

              overlaps = False
              if direction == 'right' and test_overlap(head, candidate_read, curr_contig[-1], direction, farthest_support, criteria, relaxed = km):
                # option: relaxed = km
                overlaps = True
              if direction == 'left' and test_overlap(head, curr_contig[0], candidate_read, direction, farthest_support, criteria, relaxed = km):
                overlaps = True
              if overlaps:
                good_candidates.append(head)
                criteria[head] = len(creads[head]) / 2 - 1    # Criteria = # of neighbors in 1-deg nhood
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
            # find_genomic_position(rr[hr.index(gc)], hr, rr)       # testing

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
          # filtered_good_candidates.sort(key = lambda d: criteria[d])
          # random.shuffle(filtered_good_candidates)    # testing
          print filtered_good_candidates

          print 'TIME for candidates:', datetime.datetime.now() - curr_time
          # Once we choose a particular candidate
          consensus_temp = ''
          for i in range(len(filtered_good_candidates)):
            best_head = filtered_good_candidates[i]
            h = best_head
            print 'CANDIDATE CHOSEN:'
            test_overlap(h, rr[hr.index(h)], curr_contig[-1], direction, farthest_support, criteria, relaxed = False, print_alignment = True)     # testing
            # find_genomic_position(rr[hr.index(h)], hr, rr)  # testing
            if use_ecs and h in ecs:
              consensus_temp = ecs[h]
            else:
              consensus_temp = error_correct(ec_tool, h, headers, creads, hr, rr)
            if len(consensus_temp) != 0 and consensus_temp not in curr_contig:
              break  
          if len(consensus_temp) == 0:
            print 'COULD NOT ERROR CORRECT ANY FILTERED GOOD CANDIDATES'
            h = old_h
          else:
            print 'New header:', h, criteria[h]       # testing
            print 'SUCCESS!',                         # testing 
            con_pos = find_genomic_position(consensus_temp, hr, rr, align_consensus = True)     # testing
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


def find_genomic_position(read, hr, rr, print_alignment = False, align_consensus = False):
  if read in rr:
    head = hr[rr.index(read)]
  else:
    head = '>none'
  temp_file = 'temp_read' + temp_sig + '.fasta'
  with open(temp_file, 'w') as f:
    f.write(head + '\n' + read)

  new_blasr_options = blasr_options
  temp_blasr_options = ''
  if print_alignment:
    temp_blasr_options = '-bestn 1 -m 0'
  if align_consensus:
    temp_blasr_options += ' -maxMatch 20'
    new_blasr_options += ' -maxMatch 20'

  if print_alignment:
    print commands.getstatusoutput(blasr_exe + ' ' + temp_file +' ' + e_coli_genome + ' ' + temp_blasr_options)[1]
  status = commands.getstatusoutput(blasr_exe + ' ' + temp_file +' ' + e_coli_genome + ' ' + new_blasr_options)[1]

  if len(status.split()) > blasr_zero_len:
    print status
    acc = float(status.split()[blasr_zero + 5])
    beg = int(status.split()[blasr_zero + 6])
    end = int(status.split()[blasr_zero + 7])
    print status
    print '\taligned to:', beg, end, acc
    return (beg + end) / 2
  else:
    print '\tFAILED ALIGNMENT'
    return -1
  

def test_overlap(head1, seq1, seq2, direction, farthest_support, criteria, relaxed = False, print_alignment = False):
  # Tests that seq1 is after seq2
  # farthest_support is a list that will contains distances 
  # from the end (depending on direction) of the current read

  dist_from_end = 2000
  # dist_from_end = 100
  acc_cutoff = overlap_accuracy_cutoff
  len_cutoff = overlap_length_cutoff

  temps1 = 'temp_seq1' + temp_sig + '.fasta'
  temps2 = 'temp_seq2' + temp_sig + '.fasta'
  with open(temps1, 'w') as f:
    f.write('>1\n' + seq1)
  with open(temps2, 'w') as f:
    f.write('>2\n' + seq2)

  if print_alignment:
    temp_blasr_options = '-bestn 1 -m 1'
    print commands.getstatusoutput(blasr_exe + ' ' + temps1 + ' ' + temps2 + ' ' + temp_blasr_options)[1]

  status = commands.getstatusoutput(blasr_exe + ' ' + temps1 + ' ' + temps2 + ' ' + blasr_options)[1]
  if len(status.split()) == blasr_zero_len:
    return False
  # print status                    # TESTING
  r2_strand_dir = int(status.split()[blasr_zero + 2])
  r1_strand_dir = int(status.split()[blasr_zero + 3])
  accuracy = float(status.split()[blasr_zero + 5])
  beg_align_r1 = int(status.split()[blasr_zero + 6])
  end_align_r1 = int(status.split()[blasr_zero + 7])
  total_len_r1 = int(status.split()[blasr_zero + 8])
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[blasr_zero + 9])
  end_align_r2 = int(status.split()[blasr_zero + 10])
  total_len_r2 = int(status.split()[blasr_zero + 11])
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2          # Average alignment length

  # update criteria
  # criteria[head1] = length
  criteria[head1] = accuracy

  if r2_strand_dir != r1_strand_dir:
    return False

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

  # Traditional 1-deg nhood
  for k in ktmers:
    next_read = find_extending_read(k, headers, hr, rr)
    if len(next_read) != 0:
      accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
      reads += accepted

  # Bounded Nhood
  # reads = [s for s in get_nhood(header, headers, creads) if s not in traversed_headers]

  if len(reads) != 0:
    return reads

  return reads

  # # Try bounded n-deg nhood
  # print 'trying n-deg bounded nhood'
  # reads = [s for s in get_nhood(header, headers, creads) if s not in traversed_headers]
  # print 'found', len(reads), 'possible reads'
  # return reads

  # # Try finding a read that comes close to passing the current read
  # for k in dist_to_end.keys():
  #   dist_to_end[k] -= leniency

  # print 'trying n-degree nhood extension'
  # traversed_ktmers = set(ktmers)
  # while True:
  #   print '\tkmers considered:', len(traversed_ktmers)
  #   new_ktmers = defaultdict(list)   # Key = new-ktmer, Val = [old ktmer that is connected]
  #   if len(ktmers) > num_kmers_cutoff:
  #     print 'Stopping - too many kt-mers'
  #     return ''
  #   for kt in ktmers:
  #     for kn in find_neighboring_ktmers(kt, headers, creads):
  #       if kn not in traversed_ktmers:
  #         if kn not in new_ktmers.keys() or kt not in new_ktmers[kn]:
  #           new_ktmers[kn].append(kt)
  #   if len(new_ktmers.keys()) == 0:
  #     print 'None found'
  #     return ''

  #   print 'found', len(new_ktmers.keys()), 'ktmers from', len(ktmers), 'ktmers'
  #   new_dist_to_end = dict()
  #   reads = []
  #   num_neighbors = []
  #   for kt in new_ktmers.keys():
  #     k = new_ktmers[kt][0]
  #     dist = dist_bw_ktmers(kt, k, headers, creads)
  #     if dist == None:
  #       continue
  #     if direction == 'right':
  #       new_dist_to_end[kt] = dist_to_end[k] + dist
  #     if direction == 'left':
  #       new_dist_to_end[kt] = dist_to_end[k] - dist
  #     next_read = find_extending_read(k, headers, hr, rr)
  #     if len(next_read) != 0:
  #       accepted = [nr for nr in next_read if nr not in traversed_headers and nr not in reads]
  #       reads += accepted
  #       for nr in next_read:
  #         num_neighbors.append(len(creads[nr]) / 2 - 1)
  #   if len(reads) != 0:
  #     return [x for (y, x) in sorted(zip(num_neighbors, reads), reverse = True)]

  #   # Filter out those who are too far
  #   for k in ktmers:
  #     traversed_ktmers.add(k)

  #   # print new_dist_to_end
  #   ktmers = []
  #   for kt in new_ktmers.keys():
  #     if new_dist_to_end[kt] < backtrack_limit and kt not in traversed_ktmers:
  #       ktmers.append(kt)

  #   dist_to_end = new_dist_to_end


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


def get_nhood(header, headers, creads):
  def find_acceptable_neighbors(ktmer, header, creads, left, right):
    ktmers = []
    dists = dict()
    kt_pos = -1
    temp_total = 0
    for i in range(len(creads[header])):
      if i % 2 == 1:
        ktmers.append(creads[header][i])
        dists[creads[header][i]] = temp_total
        if creads[header][i] == ktmer:
          kt_pos = temp_total
      else:
        temp_total += int(creads[header][i])

    for key in dists.keys():
      dists[key] -= kt_pos
      if -left <= dists[key] <= right and key != ktmer:
        ktmers.append(key)
    return ktmers, dists

  ktmers = []
  nhood_headers = []
  dist_to_left = dict()     # Key = ktmer, Val = int dist
  dist_to_right = dict()     # Key = ktmer, Val = int dist
  temp_total = 0

  for i in range(len(creads[header])):
    if i % 2 == 1:
      kt = creads[header][i]
      ktmers.append(kt)
      for h in headers[kt]:
        if h not in nhood_headers:
          nhood_headers.append(h)
      dist_to_left[kt] = temp_total
    else:
      temp_total += int(creads[header][i])
  for kt in dist_to_left.keys():
    dist_to_right[kt] = temp_total - dist_to_left[kt]

  # print len(nhood_headers)
  num_iterations = 0
  used_ktmers = set(ktmers)
  exit = False
  print 'starting w/', len(nhood_headers), 'headers'
  while len(nhood_headers) < nhood_header_limit and num_iterations < nhood_it_limit:
    num_iterations += 1
    print num_iterations
    # print 'loop'              # testing
    # new_ktmers = copy.copy(ktmers)
    new_ktmers = []
    for kt in ktmers:
      for h in headers[kt]:
        kts_new, dists = find_acceptable_neighbors(kt, h, creads, dist_to_left[kt], dist_to_right[kt])
        for kt_new in kts_new:
          if kt_new not in ktmers and kt_new not in used_ktmers:
            new_ktmers.append(kt_new)
            dist_to_left[kt_new] = dist_to_left[kt] + dists[kt_new]
            dist_to_right[kt_new] = dist_to_right[kt] - dists[kt_new]
            nhood_headers += [s for s in headers[kt_new] if s not in nhood_headers]
            if len(nhood_headers) >= nhood_header_limit:
              exit = True
              break
        if exit:
          break
      if exit:
        break
      print len(nhood_headers), len(new_ktmers), len(used_ktmers)              # testing
    for s in new_ktmers:
      used_ktmers.add(s)
    ktmers = copy.copy(new_ktmers)

  return nhood_headers

def get_1_deg_nhood(header, creads, headers):
  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    print header
    return ''
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      if h != header:
        collected_h.add(h)
      # find_genomic_position(rr[hr.index(h)], hr, rr)    # testing
  return collected_h

def error_correct(ec_tool, header, headers, creads, hr, rr, temp_sig_out = None):
  if temp_sig_out is not None:
    temp_sig = temp_sig_out
  else:
    global temp_sig
  reads = []

  print 'POSITION OF BASE READ'                     # testing
  find_genomic_position(rr[hr.index(header)], hr, rr)       # testing
  
  # 1-deg nhood
  collected_h = get_1_deg_nhood(header, creads, headers)

  # print 'FINDING POSITIONS OF 1-DEG NHOOD READS'  # testing
  # for ch in collected_h:                          # testing
    # find_genomic_position(rr[hr.index(ch)], hr, rr)       # testing

  # n-degree nhood
  # collected_h = get_nhood(header, headers, creads)

  # Use 2-deg nhood, no width bound (irrelevant reads, but ec tool should handle)
  # new_ktmers = []
  # new_collected_h = copy.copy(collected_h)
  # for ch in collected_h:
  #   for i in range(len(creads[ch])):
  #     if i % 2 == 1:
  #       new_ktmers.append(creads[ch][i])
  #   for kt in new_ktmers:
  #     for h in headers[kt]:
  #       new_collected_h.add(h)
  # collected_h = new_collected_h

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
    text = f.readlines()
    if len(text) > 1:
      consensus = text[1].strip()
      if len(consensus) == 0:
        return ''
    else:
      consensus = ''
  print 'consensus len:', len(consensus), 'out of', len(rr[hr.index(header)])
  return consensus.upper()


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
  for i in range(len(h)):
    h[i] = h[i].split()[0]
  with open(creads_file) as f:
    for i, line in enumerate(f):
      creads[h[i]] = line.split()
  return creads


def build_headers_dict(ktmer_headers_file):
  headers = defaultdict(list)   # Key = ktmer, Val = [headers]
  unique = set()
  with open(ktmer_headers_file) as f:
    for i, line in enumerate(f):
      words = line.split()
      headers[words[0]] = words[1:]
      for w in words[1:]:
        unique.add(w)
  print 'Found', len(unique), 'reads with kt-mers'
  return headers



def find_jumps_in_contigs(contigs_fold, parallel_prefix):
  sc_overlap_len = 10000
  sc_overlap_acc = 90
  dist_to_end = 200
  jump_pct_cutoff = 0.60        # If > this pct of candidates show jumps, then split up contig
  min_support = 3               # When splitting contigs, split at any location with this much support

  # Align one contig against all others
  traversed = set()
  combined_files = [s for s in os.listdir(contigs_fold) if fnmatch.fnmatch(s, '*combined.fasta')]
  combined_files = [s for s in combined_files if s not in traversed]
  for fn in combined_files:
    if fn.split('_')[1][:1] == parallel_prefix:
      curr_contig = fn
      cc_file = contigs_fold + curr_contig
      traversed.add(fn)
    else:
      continue

    print cc_file

    jump_pos = dict()   # Key = genomic position (int), Val = count
    jump_denom = 0
    jump_score = 0
    for fn2 in os.listdir(contigs_fold):
      if fnmatch.fnmatch(fn2, '*combined.fasta') and fn2 != curr_contig:
        fn_file = contigs_fold + fn2
        status = commands.getstatusoutput(blasr_exe + ' ' + cc_file + ' ' + fn_file + ' ' + blasr_options)[1]
        if len(status.split()) == blasr_zero_len:
          # print 'FAILED BLASR ALIGNMENT'    # testing
          continue
        else:
          # print status                      # TESTING
          acc = float(status.split()[blasr_zero + 5])
          beg_align_r1 = int(status.split()[blasr_zero + 6])
          end_align_r1 = int(status.split()[blasr_zero + 7])
          total_len_r1 = int(status.split()[blasr_zero + 8])
          end_pos_r1 = total_len_r1 - end_align_r1
          beg_align_r2 = int(status.split()[blasr_zero + 9])
          end_align_r2 = int(status.split()[blasr_zero + 10])
          total_len_r2 = int(status.split()[blasr_zero + 11])
          end_pos_r2 = total_len_r2 - end_align_r2
          length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2

          if length > sc_overlap_len and acc > sc_overlap_acc:
            jump_denom += 1
            print status
            if beg_align_r1 > dist_to_end and beg_align_r2 > dist_to_end or end_pos_r1 > dist_to_end and  end_pos_r2 > dist_to_end:
              print 'found jump'
              jump_score += 1
              if beg_align_r2 not in jump_pos:
                jump_pos[beg_align_r2] = 0
              if end_align_r2 not in jump_pos:
                jump_pos[end_align_r2] = 0
              jump_pos[beg_align_r2] += 1
              jump_pos[end_align_r2] += 1
    print 'jump score:', jump_pos, jump_denom

    # Smooth out counts on either end of contig, from 0 to dist_to_end, and dist_to_end to last bp
    total_contig_len = total_len_r2
    new_jump_pos = {0: 0, total_contig_len: 0}
    for k in jump_pos:
      if k < dist_to_end:
        new_jump_pos[0] += jump_pos[k]
      if k > total_contig_len - dist_to_end:
        new_jump_pos[total_contig_len] += jump_pos[k]
      if dist_to_end <= k <= total_contig_len - dist_to_end:
        if k not in new_jump_pos:
          new_jump_pos[k] = jump_pos[k]
        else:
          new_jump_pos[k] += jump_pos[k] 
    jump_pos = new_jump_pos
    print 'after smoothing:', jump_pos, jump_denom

    # Split
    if jump_denom > 0 and float(jump_score) / float(jump_denom) >= jump_pct_cutoff:
      split_pts = []
      for k in jump_pos:
        if jump_pos[k] >= min_support:
          split_pts.append(k)
      split_pts.sort()
      if split_pts[0] != 0:
        split_pts.insert(0, 0)
      if split_pts[-1] != total_contig_len:
        split_pts.append(total_contig_len)
    else:
      split_pts = [0, total_contig_len]

    # Write split contigs
    ch, cr = rf.read_fasta(contigs_fold + curr_contig)
    ch = ch[0]
    cr = cr[0]
    for i in range(len(split_pts) - 1):
      # [:-6] is to remove '.fasta'
      new_contig_file = contigs_fold + curr_contig[: -6] + '_split_' + str(i) + '.fasta'
      with open(new_contig_file, 'w') as f:
        f.write(ch + '_split_' + str(i) + '\n' + cr[split_pts[i] : split_pts[i + 1]])




def output_all_1_deg_nhoods(reads_file, creads_file, ktmer_headers_file, ec_tool, parallel_prefix):
  out_fold = '/home/mshen/research/1deg_nhoods_20kb_c30full/'
  hr, rr = rf.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)

  par_range = range(len(hr)) 
  if parallel_prefix == '000':
    par_range = range(3000)
  if parallel_prefix == '100':
    par_range = range(3000, 6000)
  if parallel_prefix == '200':
    par_range = range(6000, 9000)
  if parallel_prefix == '300':
    par_range = range(9000, 12000)
  if parallel_prefix == '400':
    par_range = range(12000, 15000)
  if parallel_prefix == '500':
    par_range = range(15000, 18000)
  if parallel_prefix == '600':
    par_range = range(18000, len(hr))

    
  for i in par_range:
    print i
    header = hr[i]
    collected_h = get_1_deg_nhood(header, creads, headers)

    base_file = str(i) + '_base.fasta'
    hood_file = str(i) + '_hood.fasta'
    with open(base_file, 'w') as f:
      f.write(header + '\n' + rr[i])
    with open(hood_file, 'w') as f:
      for h in collected_h:
        f.write(h + '\n' + rr[hr.index(h)] + '\n')

    ec_out = ec_prefix + base_file
    new_ec_out = str(i) + '_corr.fasta'
    status = commands.getstatusoutput(ec_tool + ' ' + base_file + ' ' + hood_file)[1]
    commands.getstatusoutput('mv ' + ec_out + ' ' + new_ec_out)
    commands.getstatusoutput('mv ' + new_ec_out + ' ' + out_fold)
    commands.getstatusoutput('mv ' + base_file + ' ' + out_fold)
    commands.getstatusoutput('mv ' + hood_file + ' ' + out_fold)
  return


def output_some_1_deg_nhoods(contigs_results_file, reads_file, creads_file, ktmer_headers_file, ec_tool):
  out_fold = '/home/mshen/research/1deg_nhoods_20kb_c30full/'
  hr, rr = rf.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  
  input_headers = []
  with open(contigs_results_file) as f:
    for i, line in enumerate(f):
      h = '/'.join(line.split()[0].split('/')[:-1])
      if h[-5:] == 'START':
        h = h[:-5]
      print h
      input_headers.append(h)

  creads = build_creads_dict(creads_file, reads_file)
  headers = build_headers_dict(ktmer_headers_file)



  for i in range(len(input_headers)):
    header = input_headers[i]
    collected_h = get_1_deg_nhood(header, creads, headers)

    base_file = str(i) + '_base.fasta'
    hood_file = str(i) + '_hood.fasta'
    with open(base_file, 'w') as f:
      f.write(header + '\n' + rr[i])
    with open(hood_file, 'w') as f:
      for h in collected_h:
        f.write(h + '\n' + rr[hr.index(h)] + '\n')

    ec_out = ec_prefix + base_file
    new_ec_out = str(i) + '_corr.fasta'
    status = commands.getstatusoutput(ec_tool + ' ' + base_file + ' ' + hood_file)[1]
    commands.getstatusoutput('mv ' + ec_out + ' ' + new_ec_out)
    commands.getstatusoutput('mv ' + new_ec_out + ' ' + out_fold)
    commands.getstatusoutput('mv ' + base_file + ' ' + out_fold)
    commands.getstatusoutput('mv ' + hood_file + ' ' + out_fold)
  return


def read_ec_from_file():
  # Takes 13 seconds to read 24k files
  ec_fold = '/home/mshen/research/1deg_nhoods/'
  ecs = dict()      # Key = header, Val = consensus
  
  print 'Building EC consensus dictionary from', ec_fold, '...'
  for fn in os.listdir(ec_fold):
    if fnmatch.fnmatch(fn, '*corr.fasta'):
      with open(ec_fold + fn) as f:
        lines = f.readlines()
      ecs[lines[0].strip()] = lines[1].strip()
  return ecs


def check_contigs(contigs_fold, reads_file):
  # Aligns component reads to combined contigs in an effort to find jumps.
  # 1/5/15: Doesn't work very well.
  hr, rr = rf.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
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
  dist_to_end = 2000
  blasr_options = '-bestn 1 -m 1 -maxMatch 20'
  temp_base = 'temp_contigbase_' + temp_sig + '.fasta'
  temp_try = 'temp_contigtry_' + temp_sig + '.fasta'
  for fn in os.listdir(contigs_fold):
    if fnmatch.fnmatch(fn, '*[0-9].fasta'):
      print fn                              # TESTING
      # if fn.split('.')[0] + '_combined.fasta' in os.listdir(contigs_fold):
        # continue
      hs, rs = rf.read_fasta(contigs_fold + fn)
      for i in range(len(hs)):
        hs[i] = hs[i].split()[0]
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
          if len(status.split()) == blasr_zero_len:
            print 'FAILED BLASR ALIGNMENT'
            num_fails += 1
            new_bases.append(new_base)
            continue
          else:
            print status, '\n', len(status.split())                      # TESTING
            acc = float(status.split()[blasr_zero + 5])
            beg_align_r1 = int(status.split()[blasr_zero + 6])
            end_align_r1 = int(status.split()[blasr_zero + 7])
            total_len_r1 = int(status.split()[blasr_zero + 8])
            end_pos_r1 = total_len_r1 - end_align_r1
            beg_align_r2 = int(status.split()[blasr_zero + 9])
            end_align_r2 = int(status.split()[blasr_zero + 10])
            total_len_r2 = int(status.split()[blasr_zero + 11])
            end_pos_r2 = total_len_r2 - end_align_r2

            # If we wrap around
            # START TERRIBLE COPY/PASTE
            if beg_align_r1 == 0 and acc > acc_cutoff:
              with open(temp_base, 'w') as f:
                f.write('>base\n' + base[-10000:])
              with open(temp_try, 'w') as f:
                f.write('>try\n' + rs[i])
              status = commands.getstatusoutput(blasr_exe + ' ' + temp_try + ' ' + temp_base + ' ' + blasr_options)[1]
              if len(status.split()) == blasr_zero_len:
                print 'FAILED BLASR ALIGNMENT'
                num_fails += 1
                new_bases.append(new_base)
                continue
              else:
                print status                      # TESTING
                acc = float(status.split()[blasr_zero + 5])
                beg_align_r1 = int(status.split()[blasr_zero + 6])
                end_align_r1 = int(status.split()[blasr_zero + 7])
                total_len_r1 = int(status.split()[blasr_zero + 8])
                end_pos_r1 = total_len_r1 - end_align_r1
                beg_align_r2 = int(status.split()[blasr_zero + 9])
                end_align_r2 = int(status.split()[blasr_zero + 10])
                total_len_r2 = int(status.split()[blasr_zero + 11])
                end_pos_r2 = total_len_r2 - end_align_r2

              if end_pos_r1 < dist_to_end: # and beg_align_r2 < dist_to_end:
                if acc > acc_cutoff:
                  new_base = base[: len(base) - end_pos_r1]
                  new_base += rs[i][end_align_r2 : ]

            # END TERRIBLE COPY/PASTE

            elif end_pos_r1 < dist_to_end: # and beg_align_r2 < dist_to_end:
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
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]

  for kt in headers.keys():
    print kt
    clusters = []
    for h in headers[kt]:
      tempfile = 'temp' + temp_sig + '.fasta'
      with open(tempfile, 'w') as f:
        f.write('>1\n' + rr[hr.index(h)]) 
      status = commands.getstatusoutput(blasr_exe + ' ' + tempfile + ' ' + e_coli_genome + ' ' + blasr_options)[1]
      if len(status.split()) > blasr_zero_len:
        beg = int(status.split()[blasr_zero + 6])
        end = int(status.split()[blasr_zero + 7])
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
    if len(status.split()) != blasr_zero_len:
      acc = float(status.split()[blasr_zero + 5])
      beg_align_r1 = int(status.split()[blasr_zero + 6])
      end_align_r1 = int(status.split()[blasr_zero + 7])
      total_len_r1 = int(status.split()[blasr_zero + 8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[blasr_zero + 9])
      end_align_r2 = int(status.split()[blasr_zero + 10])
      total_len_r2 = int(status.split()[blasr_zero + 11])
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
    if len(status.split()) != blasr_zero_len:
      acc = float(status.split()[blasr_zero + 5])
      beg_align_r1 = int(status.split()[blasr_zero + 6])
      end_align_r1 = int(status.split()[blasr_zero + 7])
      total_len_r1 = int(status.split()[blasr_zero + 8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[blasr_zero + 9])
      end_align_r2 = int(status.split()[blasr_zero + 10])
      total_len_r2 = int(status.split()[blasr_zero + 11])
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


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
