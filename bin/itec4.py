# Treats reads as the basis for error correction and extension, rather than kt-mers as in iterative_ec3.py
#
# Related files:
# combine_ec_contigs.py : Combines the read sets produced by this code into contigs

import sys, string, datetime, random, copy, os, commands, fnmatch
import numpy as np
from collections import defaultdict

import mylib as ml
import find_read
import kmer_matching
import convert_creads_to_nhoods
import trim_circular

global temp_sig
global contigs_fold
temp_sig = str(datetime.datetime.now()).split()[1]

PRIOR = '/home/yu/max/research/'
OVERLAP_ACCURACY_CUTOFF = 75    # .
OVERLAP_LENGTH_CUTOFF = 7000     # .
OVERLAP_ACCURACY_CUTOFF_CONSENSUS = 98
OVERLAP_LENGTH_CUTOFF_CONSENSUS = 7000
MIN_EXTENSION = 1               # Min. bp extension candidates need to extend
# OVERLAP_LENGTH_CUTOFF = 300     # .
NUM_ATTEMPTS = 1                # Number of times to try nhood extension.
SUPPORT_CUTOFF = 70             # CANDIDATE: Required pct accuracy for support to count
SUPPORT_RATIO = 0.6             # CANDIDATE: Required support for a chosen read from other candidates
LIMIT_KM_TIMES_TOTAL = 3        # How many times to attempt k-mer matching extension per direction
KM_K = 15                       # .
KM_CUTOFF = 100                 # .
SUPPORT_DIST_CUTOFF = 100000    # CONSENSUS: Bp. length, acceptable support distance from end of consensus
SUPPORT_T = 3                   # CONSENSUS: Req. # reads to support a position to determine farthest support
NHOOD_HEADER_LIMIT = float('inf')         # .
N21RATIO_CUTOFF = 0.05             # If n2/n1 is greater than this, find another consensus
# BLASR_EXE = '/home/jeyuan/blasr/alignment/bin/blasr'
BLASR_EXE = 'blasr'
BLASR_ZERO = 4      # 0 on debruijn, 4 on Yu's computer
BLASR_ZERO_LEN = 8  # 0 on debruijn, 8 on Yu's computer
BLASR_OPTIONS = '-bestn 1 -m 1'   # Concise output
# E_COLI_GENOME = '/home/mshen/research/data/E_COLI_GENOME.fasta'
# E_COLI_GENOME = '/home/yu/E_COLI_GENOME.fasta'
E_COLI_GENOME = '/home/yu/data/ecoli_consensus_mark.fasta'
USE_ECS = False

# IMPORTANT - CHANGE THIS WHEN CHANGING EC TOOL
# EC_PREFIX = 'C0421_'
EC_PREFIX = 'C0706_'
EC_N = '10'  # CHANGE THIS DEPENDING ON COVERAGE, YU RECOMMENDS = 0.2 * COV

def main():
  global contigs_fold
 
  cf_dir = sys.argv[1]
  contigs_fold = PRIOR + cf_dir
  parallel_prefix = sys.argv[6]
  # parallel_prefix = '0'
  cov = sys.argv[2]
  _k = sys.argv[3]
  _t = sys.argv[4]
  _num = sys.argv[5]

  if not os.path.exists(contigs_fold):
    os.makedirs(contigs_fold)

  # NEW 20KB DATASET
  # reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  # reads_file = '/home/mshen/research/data/reads.20k.' + cov + 'x.rc.fasta'
  reads_file = PRIOR + 'data/undersampled_20k/reads.20k.' + cov + 'x.' + _num + 'rc.fasta'
  # reads_file = PRIOR + 'data/reads.20k.rc.fasta'

  # creads_file = '/home/mshen/research/data/temp_creads.out' + cov + 'x_' + _k + '_' + _t + '_rc.out'
  # ktmer_headers_file = '/home/mshen/research/data/temp_ktmer_headers' + cov + 'x_' + _k + '_' + _t + '_rc.out'
  creads_file = PRIOR + 'data/undersampled_20k/temp_creads.out' + cov + 'x_' + _k + '_' + _t + '_rc_v2_' + _num + '.out'
  ktmer_headers_file = PRIOR + 'data/undersampled_20k/temp_ktmer_headers' + cov + 'x_' + _k + '_' + _t + '_rc_v2_' + _num + '.out'

  # creads_file = PRIOR + 'data/temp_creads.outrx_27_6_rc_v2.out'
  # ktmer_headers_file = PRIOR + 'data/temp_ktmer_headersrx_27_6_rc_v2.out'

  # ec_tool = '/home/yu/program/error_correction_0701.sh'
  ec_tool = '/home/yu/program/error_correction_0706.sh'
  # ec_tool = '/home/yu/program/' + sys.argv[2]

  print 'Reads File:', reads_file, '\ncreads File:', creads_file, '\nktmer Headers File:', \
    ktmer_headers_file, '\nEC Tool:', ec_tool, '\nContigs fold', contigs_fold
  print 'Cov, k, t, num', cov, _k, _t, _num


  # Actions
  iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix)
  # combine_contigs(contigs_fold)
  # ktmer_reads_pct_overlap(ktmer_headers_file, reads_file)
  # check_contigs(contigs_fold, reads_file)
  # output_all_1_deg_nhoods(reads_file, creads_file, ktmer_headers_file, ec_tool, parallel_prefix)
  # contigs_results_file = '/home/mshen/research/contigs30/contig_70results.fasta'
  # output_some_1_deg_nhoods(contigs_results_file, reads_file, creads_file, ktmer_headers_file, ec_tool)
  # find_jumps_in_contigs(contigs_fold, parallel_prefix)
  # convert_creads_to_nhoods.convert_creads_to_nhoods(reads_file, creads_file, ktmer_headers_file)


def iterative_ec(reads_file, ktmer_headers_file, creads_file, ec_tool, parallel_prefix):
  global contigs_fold
  hr, rr = ml.read_fasta(reads_file)
  headers = build_headers_dict(ktmer_headers_file)
  creads = build_creads_dict(creads_file, hr, rr)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  ktmers = headers.keys()
  random.shuffle(ktmers)
  print 'Found', len(ktmers), 'kt-mers.'
  # ktmers = filter_ktmers(ktmers, creads, headers)
  # print 'After filtering,', len(ktmers), 'kt-mers remain.'

  completed_contigs = []    # Holds combined contigs
  past_headers = set()
  num_contig_attempts = 2                   # testing
  print 'Trying to make', num_contig_attempts, 'contigs'
  # num_contig_attempts = len(ktmers)
  for m in range(num_contig_attempts):
    print '\n' + str(datetime.datetime.now())
    h = find_new_read(ktmers, completed_contigs, past_headers, headers, creads)
    print 'STARTING HEADER:\n', h

    scon = error_correct(ec_tool, h, headers, creads, hr, rr)[0]
    curr_contig = [scon]
    ccc = scon
    if len(scon) == 0:
      print 'Empty consensus'
      continue
    if find_in_contigs(completed_contigs, scon):
      print 'Found consensus in previous contig'
      continue 
    
    
    curr_contig_headers = [h]
    curr_contig_headers_data = ['_']
    master_h = h
    master_traversed_headers = [h]
    curr_time = datetime.datetime.now()

    # MAIN LOOP
    # for direction in ['right', 'left']:
    for direction in ['right']:
      counter = 0
      limit_km_times = LIMIT_KM_TIMES_TOTAL
      h = master_h
      while True:
        # Break condition: Current header doesn't change, meaning we couldn't find any extension candidates
        counter += 1
        print datetime.datetime.now() - curr_time, datetime.datetime.now()
        print '-----------------\niteration', counter, direction
        curr_time = datetime.datetime.now()
        if counter > 1500:
          print 'Reached maximum number of iterations:', counter
          break
        old_h = h
        if limit_km_times > 0:
          NUM_ATTEMPTS_temp = NUM_ATTEMPTS + 1
        else:
          NUM_ATTEMPTS_temp = NUM_ATTEMPTS
        for i in range(NUM_ATTEMPTS_temp):
          print 'Attempt', i,                                            # testing
          km = False
          km_early_out = False
          traversed_headers = master_traversed_headers

          # Grab candidates via nhood extension or kmer matching
          if i < NUM_ATTEMPTS:
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
              possible_heads = kmer_matching.kmer_matching(curr_contig[-1], reads_file, KM_K, KM_CUTOFF, file_bool = False)
            if direction == 'left':
              possible_heads = kmer_matching.kmer_matching(curr_contig[0], reads_file, KM_K, KM_CUTOFF, file_bool = False)
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
            candreads = []
            for head in possible_heads:
              candreads.append(rr[hr.index(head)])

            if direction == 'right':
              good_candidates = test_overlap_multiple(possible_heads, candreads, curr_contig[-1], direction, farthest_support, criteria, relaxed = km)
            else:
              good_candidates = test_overlap_multiple(possible_heads, candreads, curr_contig[0], direction, farthest_support, criteria, relaxed = km)

            for gc in good_candidates:
              criteria[gc] = len(creads[gc]) / 2 - 1

            if len(farthest_support) == 0:
              break
            support_pos = min(SUPPORT_T - 1, len(farthest_support) - 1)
            farthest_support.sort()
            best_support = farthest_support[support_pos]
            if best_support < SUPPORT_DIST_CUTOFF:
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

          # Filter candidates that are already in curr_contig_headers
          filtered_good_candidates = [s for s in good_candidates if s not in curr_contig_headers]
          print 'Filtered', len(good_candidates) - len(filtered_good_candidates), 'old reads'
          if len(filtered_good_candidates) == 0:
            print 'No reads passed filtering'
            if km:
              km_early_out = True
            continue

          # Sort candidates by some criteria
          filtered_good_candidates.sort(key = lambda d: criteria[d], reverse = True)
          # filtered_good_candidates.sort(key = lambda d: criteria[d])
          # random.shuffle(filtered_good_candidates)    # testing
          # print filtered_good_candidates          # testing

          print 'TIME for candidates:', datetime.datetime.now() - curr_time

          # Find consensus and extend
          consensus_temp = ''
          for i in range(len(filtered_good_candidates)):
            best_head = filtered_good_candidates[i]
            h = best_head
            print 'CANDIDATE CHOSEN:'

            # find_genomic_position(rr[hr.index(h)], hr, rr)  # testing
            consensus_temp, n1, n2 = error_correct(ec_tool, h, headers, creads, hr, rr, candidates = filtered_good_candidates)

            if len(consensus_temp) != 0:
              if direction == 'right' and test_overlap(h, consensus_temp, curr_contig[-1], direction, farthest_support, criteria, print_alignment = True, consensus = True):
                break
              if direction == 'left' and test_overlap(h, curr_contig[0], consensus_temp, direction, farthest_support, criteria, print_alignment = True, consensus = True):
                break
          if len(consensus_temp) == 0:
            print 'COULD NOT ERROR CORRECT ANY FILTERED GOOD CANDIDATES'
            h = old_h
          else:
            ccc, change = extend_attach(ccc, consensus_temp, direction)
            if change == False:
              print 'New consensus does not overlap with the end of current contig'
              h = old_h
              break  # No more attempts, end current contig
            if change == True:
              print 'New header:', h, criteria[h]       # testing
              print 'SUCCESS!',                         # testing 
              if direction == 'right':
                curr_contig.append(consensus_temp)
                curr_contig_headers.append(h)
                curr_contig_headers_data.append('_' + str(n1) + '_' + str(n2))
              if direction == 'left':
                curr_contig.insert(0, consensus_temp)
                curr_contig_headers.insert(0, h)
                curr_contig_headers_data.insert(0, '_' + str(n1) + '_' + str(n2))
              master_traversed_headers.append(h)
              print 'Curr contig len:', len(ccc)
              break   # Break from more attempts

        if counter % 10 == 0:
          if find_in_contigs(completed_contigs, consensus_temp):
            'Current contig is found to overlap with existing contig'
            break

        if h == old_h:
          print 'No new reads overlapped with current contig'
          break

    # ASSESS RESULTS
    print old_h
    contig = ''
    contig_file = contigs_fold + 'contig_' + parallel_prefix + str(m) + '.fasta'
    contig_result = contigs_fold + 'contig_' + parallel_prefix + str(m) + 'results.fasta'
    contig_cc = contigs_fold + 'contig_' + parallel_prefix + str(m) + '.combined.fasta'
    print 'Writing to files...'
    for j in range(len(curr_contig)):
      if curr_contig[j] != '':
        contig += '>' + curr_contig_headers[j] + curr_contig_headers_data[j] + '\n' + curr_contig[j] + '\n'
    with open(contig_file, 'w') as f:
      f.write(contig)
    with open(contig_cc, 'w') as f:
      f.write('>combined\n' + ccc)
    status = commands.getstatusoutput(BLASR_EXE + ' ' + contig_file +' ' + E_COLI_GENOME + ' ' + BLASR_OPTIONS + ' -maxMatch 20 > ' + contig_result)[1]
    print 'Trimming contig to make it circular...'
    trim_circular.trim_circular(contig_cc)

    for s in curr_contig_headers:
      past_headers.add(s)
    completed_contigs.append(ccc)
    # END LOOP


def find_new_read(ktmers, completed_contigs, past_headers, headers, creads):
  for i in range(len(ktmers)):
    kt = ktmers[i]
    for cc in completed_contigs:
      if kt in cc:
        continue
    print 'Skipped over', i, 'kt-mers'
    h = get_read_with_most_neighbors(kt, headers, creads)
    if h not in past_headers:
      return h


def find_genomic_position(read, hr, rr, print_alignment = False, align_consensus = False):
  if read in rr:
    head = hr[rr.index(read)]
  else:
    head = '>none'
  temp_file = 'temp_read' + temp_sig + '.fasta'
  with open(temp_file, 'w') as f:
    f.write(head + '\n' + read)

  new_BLASR_OPTIONS = BLASR_OPTIONS
  temp_BLASR_OPTIONS = ''
  if print_alignment:
    temp_BLASR_OPTIONS = '-bestn 1 -m 0'
  if align_consensus:
    temp_BLASR_OPTIONS += ' -maxMatch 20'
    new_BLASR_OPTIONS += ' -maxMatch 20'

  status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_file +' ' + E_COLI_GENOME + ' ' + new_BLASR_OPTIONS)[1]
  if print_alignment:
    print status

  if len(status.split()) > BLASR_ZERO_LEN:
    print status
    acc = float(status.split()[BLASR_ZERO + 5])
    beg = int(status.split()[BLASR_ZERO + 6])
    end = int(status.split()[BLASR_ZERO + 7])
    print status
    print '\taligned to:', beg, end, acc
    return (beg + end) / 2
  else:
    print '\tFAILED ALIGNMENT'
    return -1


def test_overlap(head1, seq1, seq2, direction, farthest_support, criteria, relaxed = False, print_alignment = False, consensus = False):
  # Tests that seq1 is after seq2
  # farthest_support is a list that will contains distances 
  # from the end (depending on direction) of the current read
  # Optional consensus boolean toggles the case where we're aligning one consensus to another

  if not consensus:
    dist_from_end = 5000  # there are reads that have long tips that don't align to the genome
    # dist_from_end = 100
    acc_cutoff = OVERLAP_ACCURACY_CUTOFF
    len_cutoff = OVERLAP_LENGTH_CUTOFF
  if consensus:
    dist_from_end = 100
    acc_cutoff = OVERLAP_ACCURACY_CUTOFF_CONSENSUS
    len_cutoff = OVERLAP_LENGTH_CUTOFF_CONSENSUS
  COMPLETE_OVERLAP_THRESH = 50

  temps1 = 'temp_seq1' + temp_sig + '.fasta'
  temps2 = 'temp_seq2' + temp_sig + '.fasta'
  with open(temps1, 'w') as f:
    f.write('>1\n' + seq1)
  with open(temps2, 'w') as f:
    f.write('>2\n' + seq2)

  status = commands.getstatusoutput(BLASR_EXE + ' ' + temps1 + ' ' + temps2 + ' ' + BLASR_OPTIONS)[1]
  if print_alignment:
    print status
  if len(status.split()) == BLASR_ZERO_LEN:
    return False
  # print status                    # TESTING
  r2_strand_dir = int(status.split()[BLASR_ZERO + 2])
  r1_strand_dir = int(status.split()[BLASR_ZERO + 3])
  accuracy = float(status.split()[BLASR_ZERO + 5])
  beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
  end_align_r1 = int(status.split()[BLASR_ZERO + 7])
  total_len_r1 = int(status.split()[BLASR_ZERO + 8])
  actual_len_r1 = end_align_r1 - beg_align_r1
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
  end_align_r2 = int(status.split()[BLASR_ZERO + 10])
  total_len_r2 = int(status.split()[BLASR_ZERO + 11])
  actual_len_r2 = end_align_r2 - end_align_r1
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2   # Average alignment length

  if r2_strand_dir != r1_strand_dir:
    return False

  # Update farthest support, the distance to the end of the consensus that has support from 1deg nhood 
  if not consensus:
    # update criteria
    # criteria[head1] = length
    criteria[head1] = accuracy
    
    if direction == 'right':
      farthest_support.append(end_pos_r1)
    if direction == 'left':
      farthest_support.append(beg_align_r2)

  if not relaxed:
    if accuracy >= acc_cutoff and length > len_cutoff:
      if direction == 'right':
        if end_pos_r1 < dist_from_end and end_pos_r2 > end_pos_r1 + MIN_EXTENSION and beg_align_r2 < dist_from_end:
          return True
        if abs(actual_len_r1 - total_len_r1) < COMPLETE_OVERLAP_THRESH:
          return True

      if direction == 'left':
        if end_pos_r1 < dist_from_end and beg_align_r1 > beg_align_r2 + MIN_EXTENSION and beg_align_r2 < dist_from_end:
          return True
        if abs(actual_len_r2 - total_len_r2) < COMPLETE_OVERLAP_THRESH:
          return True
  else:
    return accuracy >= acc_cutoff and length > len_cutoff


def test_overlap_multiple(heads, seqs, seq2, direction, farthest_support, criteria, relaxed = False, print_alignment = False, consensus = False):
  # Tests that seqs overlap in direction with seq2
  if not consensus:
    dist_from_end = 5000
    acc_cutoff = OVERLAP_ACCURACY_CUTOFF
    len_cutoff = OVERLAP_LENGTH_CUTOFF
  if consensus:
    dist_from_end = 100
    acc_cutoff = OVERLAP_ACCURACY_CUTOFF_CONSENSUS
    len_cutoff = OVERLAP_LENGTH_CUTOFF_CONSENSUS
  COMPLETE_OVERLAP_THRESH = 50


  temps1 = 'temp_seq1' + temp_sig + '.fasta'
  temps2 = 'temp_seq2' + temp_sig + '.fasta'
  with open(temps1, 'w') as f:
    for i in range(len(seqs)):
      f.write('>' + str(i) + '\n' + seqs[i] + '\n')
  with open(temps2, 'w') as f:
    f.write('>2\n' + seq2)

  overlaps = []
  status = commands.getstatusoutput(BLASR_EXE + ' ' + temps1 + ' ' + temps2 + ' ' + BLASR_OPTIONS)[1]
  if print_alignment:
    print status
  if len(status.split()) == BLASR_ZERO_LEN:
    return overlaps
  # print status                    # TESTING

  newstatus = status.splitlines()[1:-1]
  for ii in range(len(newstatus)):
    line = newstatus[ii]
    r2_strand_dir = int(line.split()[2])
    r1_strand_dir = int(line.split()[3])
    accuracy = float(line.split()[5])
    beg_align_r1 = int(line.split()[6])
    end_align_r1 = int(line.split()[7])
    total_len_r1 = int(line.split()[8])
    actual_len_r1 = end_align_r1 - beg_align_r1
    end_pos_r1 = total_len_r1 - end_align_r1
    beg_align_r2 = int(line.split()[9])
    end_align_r2 = int(line.split()[10])
    total_len_r2 = int(line.split()[11])
    actual_len_r2 = end_align_r2 - beg_align_r2
    end_pos_r2 = total_len_r2 - end_align_r2

    length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2   # Average alignment length

    currindex = int(line.split('/')[0])
    currhead = heads[currindex]

    if r2_strand_dir != r1_strand_dir:
      continue

    if not consensus:
      # update criteria
      criteria[currhead] = accuracy
      
      if direction == 'right':
        farthest_support.append(end_pos_r1)
      if direction == 'left':
        farthest_support.append(beg_align_r2)

      if not relaxed:
        if accuracy >= acc_cutoff and length > len_cutoff:
          if direction == 'right':
            if end_pos_r1 < dist_from_end and end_pos_r2 > end_pos_r1 + MIN_EXTENSION and beg_align_r2 < dist_from_end:
              overlaps.append(currhead)
            if abs(actual_len_r1 - total_len_r1) < COMPLETE_OVERLAP_THRESH:
              overlaps.append(currhead)

          if direction == 'left':
            if end_pos_r1 < dist_from_end and beg_align_r1 > beg_align_r2 + MIN_EXTENSION and beg_align_r2 < dist_from_end:
              overlaps.append(currhead)
            if abs(actual_len_r2 - total_len_r2) < COMPLETE_OVERLAP_THRESH:
              overlaps.append(currhead)
      else:
        if accuracy >= acc_cutoff and length > len_cutoff:
          overlaps.append(currhead)    
  return overlaps

def find_in_contigs(completed_contigs, consensus):
  LEN_CUTOFF = 10000
  for cc in completed_contigs:
    temps1 = 'temp_seq1' + temp_sig + '.fasta'
    temps2 = 'temp_seq2' + temp_sig + '.fasta'
    with open(temps1, 'w') as f:
      f.write('>1\n' + consensus)
    with open(temps2, 'w') as f:
      f.write('>2\n' + cc)

    status = commands.getstatusoutput(BLASR_EXE + ' ' + temps1 + ' ' + temps2 + ' ' + BLASR_OPTIONS)[1]
    if len(status.split()) == BLASR_ZERO_LEN:
      return False
    accuracy = float(status.split()[BLASR_ZERO + 5])
    beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
    end_align_r1 = int(status.split()[BLASR_ZERO + 7])
    beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
    end_align_r2 = int(status.split()[BLASR_ZERO + 10])
    length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2

    if length > LEN_CUTOFF:
      print 'find in contigs:\n' + status
      return True
  return False


def extend_attach(ccc, consensus_temp, direction):
  # Extends the big contig we have right now
  dist_from_end = 100
  acc_cutoff = OVERLAP_ACCURACY_CUTOFF_CONSENSUS
  len_cutoff = OVERLAP_LENGTH_CUTOFF_CONSENSUS
  COMPLETE_OVERLAP_THRESH = 50
  blasr_o = '-bestn 1 -m 1 -maxMatch 20'

  temps1 = 'temp_seq1' + temp_sig + '.fasta'
  temps2 = 'temp_seq2' + temp_sig + '.fasta'
  with open(temps1, 'w') as f:
    f.write('>1\n' + consensus_temp)
  with open(temps2, 'w') as f:
    f.write('>2\n' + ccc)

  status = commands.getstatusoutput(BLASR_EXE + ' ' + temps1 + ' ' + temps2 + ' ' + blasr_o)[1]
  if len(status.split()) == BLASR_ZERO_LEN:
    print 'Failure: No alignment'
    print 'Files: Consensus:', temps1, 'Current contig:', temps2
    return ccc, False
  r2_strand_dir = int(status.split()[BLASR_ZERO + 2])
  r1_strand_dir = int(status.split()[BLASR_ZERO + 3])
  accuracy = float(status.split()[BLASR_ZERO + 5])
  beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
  end_align_r1 = int(status.split()[BLASR_ZERO + 7])
  total_len_r1 = int(status.split()[BLASR_ZERO + 8])
  actual_len_r1 = end_align_r1 - beg_align_r1
  end_pos_r1 = total_len_r1 - end_align_r1
  beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
  end_align_r2 = int(status.split()[BLASR_ZERO + 10])
  total_len_r2 = int(status.split()[BLASR_ZERO + 11])
  actual_len_r2 = end_align_r2 - beg_align_r2
  end_pos_r2 = total_len_r2 - end_align_r2
  length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2

  print status  # testing

  change = False
  if accuracy >= acc_cutoff and length > len_cutoff:
    if direction == 'right':
      if end_pos_r1 < dist_from_end and beg_align_r2 < dist_from_end:
        ccc = ccc[: end_align_r1] + consensus_temp[end_align_r2 :]
        change = True
      if abs(actual_len_r1 - total_len_r1) < COMPLETE_OVERLAP_THRESH:
        ccc = ccc[: end_align_r1] + consensus_temp[end_align_r2 :]
        change = True
    if direction == 'left':
      if end_pos_r1 < dist_from_end and beg_align_r2 < dist_from_end:
        ccc =  consensus_temp[: end_align_r2] + ccc[end_align_r1 :]
        change = True
      if abs(actual_len_r2 - total_len_r2) < COMPLETE_OVERLAP_THRESH:
        ccc =  consensus_temp[: end_align_r2] + ccc[end_align_r1 :]
        change = True
  if not change:
    print 'Failure:', status
    print 'Files: Consensus:', temps1, 'Current contig:', temps2
  return ccc, change


def extend_n(header, headers, creads, traversed_headers, direction, hr, rr):
  # Returns possible candidates for extension based on the 1-deg or multi-deg nhood of a raw read

  leniency = 100    # If 1-degree nhood fails, we accept a read that extends beyond border - leniency
  backtrack_limit = 1000       # Don't backtrack too far
  num_kmers_cutoff = 500       # If we start considering this many kt-mers in nhood extension, stop.
  
  dist_to_end = get_dist_to_end(header, creads, direction)
  ktmers = dist_to_end.keys()
  print len(ktmers), 'ktmers in 1-deg nhood'

  reads = []

  reads, windows = get_simple_1_deg_nhood(header, creads, headers, hr)
  reads = [s for s in reads if s not in traversed_headers]

  if len(reads) != 0:
    return reads

  return reads

  # # Try bounded n-deg nhood
  # print 'trying n-deg bounded nhood'
  # reads = [s for s in get_nhood(header, headers, creads, hr) if s not in traversed_headers]
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
  # not the current read we peml.rm 1-deg nhood extension on.
  valid = []
  # if verify_against_candidates(headers[ktmer][0], headers[ktmer][1:], hr, rr):
    # valid = headers[ktmer]
  valid = headers[ktmer]  
  return valid


def nhood_stats(base_index, nhood_indices):
  if len(nhood_indices) == 0:
    print 'empty nhood'
    return
  specificity = []
  sensitivity = []
  ground_truth = '/home/mshen/research/data/genome_nhood_7k.txt'
  rest = []
  with open(ground_truth) as f:
    for i, line in enumerate(f):
      base = line.split(':')[0]
      if base == str(base_index):
        rest = ' '.join(line.split(':')[1:]).split()
        break
  if len(rest) == 0:
    print 'no true nhood'
    return
  intersect = len(set([str(s) for s in nhood_indices]).intersection(rest))
  false_reads = [s for s in nhood_indices if str(s) not in rest]
  sensitivity.append(float(intersect) / float(len(rest)))
  specificity.append(float(intersect) / float(len(nhood_indices)))

  # print 'sensitivity:', float(sum(sensitivity)) / float(len(sensitivity))
  # print 'specificity:', float(sum(specificity)) / float(len(specificity))
  print 'false:', len(false_reads)
  print 'missing:', len(rest) - intersect
  return false_reads


def get_nhood(header, headers, creads, hr):
  def len_read(cread):
    return sum([cread[s] for s in range(len(cread)) if s % 2 == 0])

  nhood_headers, windows = get_special_1_deg_nhood(header, creads, headers, hr)
  collected_headers = [nhood_headers]   # List of lists
  collected_windows = [windows]
  collected = set(nhood_headers)

  base_index = hr.index(header)
  print header, base_index, len(list(collected))
  nhood_indices = [hr.index(s) for s in list(collected)]
  nhood_stats(base_index, nhood_indices)

  depth = 1 
  for i in range(depth):
    new_headers = []
    new_windows = []
    curr_level_headers = collected_headers[-1]
    curr_level_windows = collected_windows[-1]
    # print curr_level_windows, collected_windows

    # print len(curr_level_headers), curr_level_headers
    for j in range(len(curr_level_headers)):
      curr_head = curr_level_headers[j]
      curr_window = curr_level_windows[j]
      # print curr_window
      new_nhood_headers, new_nhood_windows = get_special_1_deg_nhood(curr_head, creads, headers, hr, curr_window)
      for k in range(len(new_nhood_headers)):
        if new_nhood_headers[k] not in collected:
          collected.add(new_nhood_headers[k])
          new_headers.append(new_nhood_headers[k])
          new_windows.append(new_nhood_windows[k])
    collected_headers.append(new_headers)
    collected_windows.append(new_windows)

    print 'degree', i + 2, len(list(collected))     # testing
    base_index = hr.index(header)
    nhood_indices = [hr.index(s) for s in list(collected)]
    false_reads = nhood_stats(base_index, nhood_indices)
    # print false_reads, [hr[s] for s in false_reads]
    # print header

  # print [hr.index(s) for s in collected]        # testing
  return list(collected)


def get_simple_1_deg_nhood(header, creads, headers, hr):
  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    print header
    return '', None
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      if h != header:
        collected_h.add(h)
      # find_genomic_position(rr[hr.index(h)], hr, rr)    # testing
  return collected_h, None


def get_special_1_deg_nhood(header, creads, headers, hr, n_range = []):
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

  def filter_special_1_deg_nhood(header, nhood_headers, creads, n_range):
    # Filters a neighborhood of reads to a master_read. Searches for overlapping kmers
    #   Span of overlapping kmers must be greater than min_bp_shared
    #   Distance b/w consecutive shared kmers must be less than max_dist



    def get_relative_dist(kmer1, kmer2, cread):
      # Returns 0 if there are no elements between kmer1 and kmer2, kmer1 must be before kmer2
      ki1 = cread.index(kmer1)
      ki2 = cread.index(kmer2)
      return sum(int(cread[s]) for s in range(ki1, ki2 + 1) if s % 2 == 0)


    leniency = 100000    # for comparing relative distances between kmers
    max_dist = 100000   # If at least one read does not have a shared kmer within this distance, False
    min_bp_shared = 3500
    extend_range = 0

    new_nhood = []
    windows = []

    master_cread = creads[header]
    for candidate in nhood_headers:
      window = [-1, -1]
      cand_cread = creads[candidate]
      master_dists = []
      cand_dists = []
      prev_cand = ''

      # Assists in n-deg nhood, only look for overlapping ktmers in range
      if n_range == []:
        master_cread_range = range(len(master_cread))
      else:
        indices = convert_pos_range_to_indices(n_range, master_cread)
        master_cread_range = range(indices[0], indices[1])

      # Search for overlapping kmers between master_cread and cand_cread
      for m_kmer in [master_cread[s] for s in master_cread_range if s % 2 == 1]:
        if m_kmer in cand_cread:
          if prev_cand == '':
            # First overlapping kmer
            window[0] = max(get_pos_in_read(m_kmer, cand_cread) - extend_range, 0)
            prev_cand = m_kmer
          else:
            m_dist = get_relative_dist(prev_cand, m_kmer, master_cread)
            c_dist = get_relative_dist(prev_cand, m_kmer, cand_cread)
            if c_dist != 0 and m_dist != 0:
              cand_dists.append(c_dist)
              master_dists.append(m_dist)
          # print m_kmer,
          prev_cand = m_kmer
      if prev_cand != '':
        window[1] = min(get_pos_in_read(prev_cand, cand_cread) + extend_range, len_read(cand_cread))

      # Filter for ktmers that are too far apart
      new_master_dists = []
      new_cand_dists = []
      rs_master = 0
      rs_cand = 0
      for i in range(len(master_dists)):
        if abs(master_dists[i] - cand_dists[i]) < leniency:
          new_master_dists.append(master_dists[i] + rs_master)
          new_cand_dists.append(cand_dists[i] + rs_cand)
          rs_master = 0
          rs_cand = 0
        else:
          rs_master += master_dists[i]
          rs_cand += cand_dists[i]
      master_dists = new_master_dists
      cand_dists = new_cand_dists

      failed = False
      for s in master_dists + cand_dists:
        if s > max_dist:
          failed = True
          continue
      if sum(master_dists) < min_bp_shared or sum(cand_dists) < min_bp_shared:
        failed = True
      if failed:
        continue

      # print hr.index(candidate), master_dists, cand_dists   # testing
      # print window      # testing
      windows.append(window)
      new_nhood.append(candidate)
    print 'Filtering nhood - prev:', len(nhood_headers), ' after:', len(new_nhood)
    return new_nhood, windows

  def remove_rc_duplicate_in_headers(headers):
    # Randomly removes both (as of 6/5/15, before just removed 1) of the duplicate reads
    # (rc and normal) in a set of headers
    # We can remove randomly because we pass into Yu's EC which uses BLASR to align
    # all headers to the base read, which will automatically flip if needed
    past = []
    prefixes = []
    for h in headers:
      trimmed_h = '_'.join(h.split('_')[1:])
      prefix = h.split('_')[0]
      if trimmed_h not in past:
        past.append(trimmed_h)
        prefixes.append(prefix)
      else:
        ind = past.index(trimmed_h)
        del past[ind]
        del prefixes[ind]
    
    ans = []
    for i in range(len(past)):
      ans.append(prefixes[i] + '_' + past[i])
    return ans

  def get_pos_in_read(kmer, cread):
    if kmer in cread:
      return sum([int(cread[s]) for s in range(cread.index(kmer)) if s % 2 == 0])
    else:
      print 'error:', kmer, 'not in', cread

  # # Start actual code for get_1_deg_nhood(...)
  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    print header
    return '', None
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      if h != header:
        collected_h.add(h)
        # find_genomic_position(rr[hr.index(h)], hr, rr)    # testing
  print 'regular nhood stats:',
  # nhood_stats(hr.index(header), [hr.index(s) for s in list(collected_h)])

  # Special 1-deg nhood
  collected_h, windows = filter_special_1_deg_nhood(header, list(collected_h), creads, n_range)
  collected_h = remove_rc_duplicate_in_headers(collected_h)
  return collected_h, windows

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

  collected_h = set()
  max_dist = 3500
  min_bp_shared = 7000
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


    # Filters out reads
    failed = False
    for d in dists:
      if d > max_dist:
        failed = True
    if sum(dists) < min_bp_shared:
      failed = True
    if failed:
      continue
    else:
      collected_h.add(k)

  collected_h = list(collected_h)
  windows = []
  for ch in collected_h:
    window = [get_pos_in_read(positions[ch][0], creads[ch]), get_pos_in_read(positions[ch][-2], creads[ch])]
    windows.append(window)
    # print window

  return collected_h, windows


def keep_duplicates_only(inp):
  # Given a list, returns a new list with 1 copy of any duplicates
  # Ex: [a, a, b] -> [a]
  d = dict()
  for item in inp:
    if item not in d:
      d[item] = 1
    else:
      d[item] += 1
  new_list = []
  for key in d.keys():
    if d[key] > 1:
      new_list.append(key)
  return new_list


def error_correct(ec_tool, header, headers, creads, hr, rr, temp_sig_out = None, candidates = []):
  if temp_sig_out is not None:
    temp_sig = temp_sig_out
  else:
    global temp_sig
  reads = []

  # print 'POSITION OF BASE READ'                     # testing
  # find_genomic_position(rr[hr.index(header)], hr, rr)       # testing
  
  # # 1-deg nhood
  collected_h, windows = get_simple_1_deg_nhood(header, creads, headers, hr)
  collected_h = list(collected_h) 
  print 'Original 1deg nhood:', len(collected_h)

  # if len(candidates) > 0:
  #   for cd in candidates:
  #     new_h, windows = get_simple_1_deg_nhood(cd, creads, headers, hr)
  #     collected_h += list(new_h)
  #   collected_h = keep_duplicates_only(collected_h)
  #   print 'Duplicates only after combine:', len(collected_h)

  # print 'FINDING POSITIONS OF 1-DEG NHOOD READS'  # testing
  # for ch in collected_h:                          # testing
    # find_genomic_position(rr[hr.index(ch)], hr, rr)       # testing

  # n-degree nhood
  # collected_h = get_nhood(header, headers, creads, hr)

  # FOR PLOTTING A-BRUIJN GRAPH OF NHOODS
  # print creads[header]
  # for ch in collected_h:
  #  print creads[ch]

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

  ec_out = EC_PREFIX + temp_orig_file
  commands.getstatusoutput('rm -rf ' + ec_out)
  # print 'Consensus out file:', ec_out
  status = commands.getstatusoutput(ec_tool + ' ' + temp_orig_file + ' ' + temp_nhood_file + ' ' + EC_N)[1]
  print status
  if 'ERROR' in status or 'No such file or directory' in status:
    print status
    return '', -1, 1

  if not os.path.isfile(ec_out):
    print 'EC did not return a consensus'
    return '', -1, 1

  ch, cr = ml.read_fasta(ec_out)
  consensus = cr[0]
  header_con = ch[0]
  if len(consensus) == 0:
      return '', -1, 1

  print header_con
  if ec_tool != '/home/yu/program/error_correction_0706.sh':
    n1 = float(header_con.split('_')[2])
    n2 = float(header_con.split('_')[3].split('(')[0])
  else:
    n1 = 1
    n2 = 0
  n21ratio = n2 / n1
  print 'consensus len:', len(consensus), 'out of', len(rr[hr.index(header)]), ' error ratio (bp):', n21ratio
  if n21ratio > N21RATIO_CUTOFF:
    return '', -1, 1
  return consensus.upper(), n1, n2


def get_read_with_most_neighbors(ktmer, headers, creads):
  best_h = ''
  best_neighbors = 0
  for h in headers[ktmer]:
    num_neighbors = len(creads[h]) / 2 - 1
    # print h, num_neighbors
    if num_neighbors > best_neighbors:
      best_neighbors = num_neighbors
      best_h = h
  return best_h


def build_creads_dict(creads_file, h, r):
  creads = defaultdict(list)   # Key = header, Val = creads 
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
        status = commands.getstatusoutput(BLASR_EXE + ' ' + cc_file + ' ' + fn_file + ' ' + BLASR_OPTIONS)[1]
        if len(status.split()) == BLASR_ZERO_LEN:
          # print 'FAILED BLASR ALIGNMENT'    # testing
          continue
        else:
          # print status                      # TESTING
          acc = float(status.split()[BLASR_ZERO + 5])
          beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
          end_align_r1 = int(status.split()[BLASR_ZERO + 7])
          total_len_r1 = int(status.split()[BLASR_ZERO + 8])
          end_pos_r1 = total_len_r1 - end_align_r1
          beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
          end_align_r2 = int(status.split()[BLASR_ZERO + 10])
          total_len_r2 = int(status.split()[BLASR_ZERO + 11])
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
    ch, cr = ml.read_fasta(contigs_fold + curr_contig)
    ch = ch[0]
    cr = cr[0]
    for i in range(len(split_pts) - 1):
      # [:-6] is to remove '.fasta'
      new_contig_file = contigs_fold + curr_contig[: -6] + '_split_' + str(i) + '.fasta'
      with open(new_contig_file, 'w') as f:
        f.write(ch + '_split_' + str(i) + '\n' + cr[split_pts[i] : split_pts[i + 1]])


def output_all_1_deg_nhoods(reads_file, creads_file, ktmer_headers_file, ec_tool, parallel_prefix):
  out_fold = '/home/yu/max/research/1deg_nhoods_27.6_55x/'
  if not os.path.exists(out_fold):
    os.makedirs(out_fold)
  hr, rr = ml.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]
  creads = build_creads_dict(creads_file, hr, rr)
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

    
  for i in [s for s in par_range if s % 2 == 1]:
    print i
    header = hr[i]
    collected_h, windows = get_simple_1_deg_nhood(header, creads, headers, hr)

    base_file = str(i) + '_base.fasta'
    hood_file = str(i) + '_hood.fasta'
    with open(base_file, 'w') as f:
      f.write(header + '\n' + rr[i])
    with open(hood_file, 'w') as f:
      for h in collected_h:
        f.write(h + '\n' + rr[hr.index(h)] + '\n')

    ec_out = EC_PREFIX + base_file
    new_ec_out = str(i) + '_corr.fasta'
    status = commands.getstatusoutput(ec_tool + ' ' + base_file + ' ' + hood_file)[1]
    commands.getstatusoutput('mv ' + ec_out + ' ' + new_ec_out)
    commands.getstatusoutput('mv ' + new_ec_out + ' ' + out_fold)
    commands.getstatusoutput('mv ' + base_file + ' ' + out_fold)
    commands.getstatusoutput('mv ' + hood_file + ' ' + out_fold)
  return


def output_some_1_deg_nhoods(contigs_results_file, reads_file, creads_file, ktmer_headers_file, ec_tool):
  out_fold = '/home/mshen/research/1deg_nhoods_20kb_c30full/'
  hr, rr = ml.read_fasta(reads_file)
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

  creads = build_creads_dict(creads_file, hr, rr)
  headers = build_headers_dict(ktmer_headers_file)



  for i in range(len(input_headers)):
    header = input_headers[i]
    collected_h, windows = get_simple_1_deg_nhood(header, creads, headers, hr)

    base_file = str(i) + '_base.fasta'
    hood_file = str(i) + '_hood.fasta'
    with open(base_file, 'w') as f:
      f.write(header + '\n' + rr[i])
    with open(hood_file, 'w') as f:
      for h in collected_h:
        f.write(h + '\n' + rr[hr.index(h)] + '\n')

    ec_out = EC_PREFIX + base_file
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
  hr, rr = ml.read_fasta(reads_file)
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
    status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_file + ' ' + comb_file + ' ' + BLASR_OPTIONS)[1]
    print status


def combine_contigs(contigs_fold):
  # 7/26/15: Might not detect if one read completely contains the current contig we are building 
  acc_cutoff = 0
  dist_to_end = 2000
  BLASR_OPTIONS = '-bestn 1 -m 1 -maxMatch 20'
  temp_base = 'temp_contigbase_' + temp_sig + '.fasta'
  temp_try = 'temp_contigtry_' + temp_sig + '.fasta'
  for fn in os.listdir(contigs_fold):
    if fnmatch.fnmatch(fn, '*[0-9].fasta'):
      print fn                              # TESTING
      # if fn.split('.')[0] + '_combined.fasta' in os.listdir(contigs_fold):
        # continue
      hs, rs = ml.read_fasta(contigs_fold + fn)
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
          status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_try + ' ' + temp_base + ' ' + BLASR_OPTIONS)[1]
          if len(status.split()) == BLASR_ZERO_LEN:
            print 'FAILED BLASR ALIGNMENT'
            num_fails += 1
            new_bases.append(new_base)
            continue
          else:
            print status                      # TESTING
            acc = float(status.split()[BLASR_ZERO + 5])
            beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
            end_align_r1 = int(status.split()[BLASR_ZERO + 7])
            total_len_r1 = int(status.split()[BLASR_ZERO + 8])
            end_pos_r1 = total_len_r1 - end_align_r1
            beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
            end_align_r2 = int(status.split()[BLASR_ZERO + 10])
            total_len_r2 = int(status.split()[BLASR_ZERO + 11])
            end_pos_r2 = total_len_r2 - end_align_r2

            # If we wrap around
            # START TERRIBLE COPY/PASTE
            if beg_align_r1 == 0 and acc > acc_cutoff:
              with open(temp_base, 'w') as f:
                f.write('>base\n' + base[-10000:])
              with open(temp_try, 'w') as f:
                f.write('>try\n' + rs[i])
              status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_try + ' ' + temp_base + ' ' + BLASR_OPTIONS)[1]
              if len(status.split()) == BLASR_ZERO_LEN:
                print 'FAILED BLASR ALIGNMENT'
                num_fails += 1
                new_bases.append(new_base)
                continue
              else:
                print status                      # TESTING
                acc = float(status.split()[BLASR_ZERO + 5])
                beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
                end_align_r1 = int(status.split()[BLASR_ZERO + 7])
                total_len_r1 = int(status.split()[BLASR_ZERO + 8])
                end_pos_r1 = total_len_r1 - end_align_r1
                beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
                end_align_r2 = int(status.split()[BLASR_ZERO + 10])
                total_len_r2 = int(status.split()[BLASR_ZERO + 11])
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
      # status = commands.getstatusoutput(BLASR_EXE + ' ' + contigs_fold + out_file + ' ' + E_COLI_GENOME + ' ' + BLASR_OPTIONS)[1]
      print 'BEST:'
      print status
      trim_circular.trim_circular(contigs_fold + out_file)
  return


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
  hr, rr = ml.read_fasta(reads_file)
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]

  for kt in headers.keys():
    print kt
    clusters = []
    for h in headers[kt]:
      tempfile = 'temp' + temp_sig + '.fasta'
      with open(tempfile, 'w') as f:
        f.write('>1\n' + rr[hr.index(h)]) 
      status = commands.getstatusoutput(BLASR_EXE + ' ' + tempfile + ' ' + E_COLI_GENOME + ' ' + BLASR_OPTIONS)[1]
      if len(status.split()) > BLASR_ZERO_LEN:
        beg = int(status.split()[BLASR_ZERO + 6])
        end = int(status.split()[BLASR_ZERO + 7])
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
    status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_file +' ' + temp_chfile + ' ' + BLASR_OPTIONS)[1]
    print status                        # TESTING
    if len(status.split()) != BLASR_ZERO_LEN:
      acc = float(status.split()[BLASR_ZERO + 5])
      beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
      end_align_r1 = int(status.split()[BLASR_ZERO + 7])
      total_len_r1 = int(status.split()[BLASR_ZERO + 8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
      end_align_r2 = int(status.split()[BLASR_ZERO + 10])
      total_len_r2 = int(status.split()[BLASR_ZERO + 11])
      end_pos_r2 = total_len_r2 - end_align_r2
      length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2     # Avg alignment length

      is_support = False
      if acc > SUPPORT_CUTOFF and length > support_len:
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
  return support_pct >= SUPPORT_RATIO


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

  print commands.getstatusoutput(BLASR_EXE + ' ' + temp_base +' ' + E_COLI_GENOME + ' ' + BLASR_OPTIONS)[1]

  support = 0
  temp_chfile = 'temp_ch_' + temp_sig + '.fasta'
  for ch in cand_r[1:]:
    with open(temp_chfile, 'w') as f:
      f.write('>' + hr[rr.index(ch)] + '\n' + ch)
    status = commands.getstatusoutput(BLASR_EXE + ' ' + temp_base +' ' + temp_chfile + ' ' + BLASR_OPTIONS)[1]
    print status                        # TESTING
    print commands.getstatusoutput(BLASR_EXE + ' ' + temp_chfile +' ' + E_COLI_GENOME + ' ' + BLASR_OPTIONS)[1]
    if len(status.split()) != BLASR_ZERO_LEN:
      acc = float(status.split()[BLASR_ZERO + 5])
      beg_align_r1 = int(status.split()[BLASR_ZERO + 6])
      end_align_r1 = int(status.split()[BLASR_ZERO + 7])
      total_len_r1 = int(status.split()[BLASR_ZERO + 8])
      end_pos_r1 = total_len_r1 - end_align_r1
      beg_align_r2 = int(status.split()[BLASR_ZERO + 9])
      end_align_r2 = int(status.split()[BLASR_ZERO + 10])
      total_len_r2 = int(status.split()[BLASR_ZERO + 11])
      end_pos_r2 = total_len_r2 - end_align_r2
      length = (end_align_r2 - beg_align_r2 + end_align_r1 - beg_align_r1) / 2     # Avg alignment length

      is_support = False
      if acc > SUPPORT_CUTOFF:
        if beg_align_r1 > dist_to_end and beg_align_r2 < dist_to_end or beg_align_r1 < dist_to_end and beg_align_r2 > dist_to_end or beg_align_r1 < dist_to_end and beg_align_r2 < dist_to_end:
          if end_pos_r1 > dist_to_end and end_pos_r2 < dist_to_end or end_pos_r1 < dist_to_end and end_pos_r2 > dist_to_end or end_pos_r1 < dist_to_end and end_pos_r2 < dist_to_end: 
            is_support = True
            support += 1
            print '+1'

  support_pct = float(support) / float(len(candidate_headers) - 1)
  print 'support found:', support_pct         # TESTING
  return support_pct >= SUPPORT_RATIO


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

  # genome = '/home/mshen/research/data/E_COLI_GENOME.fasta'
  # gh, gr = ml.read_fasta(genome)
  # for kt in ktmers:
  #   if kt in gr[0]:
  #     print gr[0].index(kt)


def ktmers_from_genome(ktmers, min_bp, max_bp):
  hg, rg = ml.read_fasta(E_COLI_GENOME)
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
