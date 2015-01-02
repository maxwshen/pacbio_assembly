# Quantifies the contigs found via itec4, given a folder containing contig files and contig results files.
# Expected file names: contig_78results.fasta

import sys, string, datetime, random, copy, os, commands, fnmatch, numpy

from collections import defaultdict

def main():
  input_folder = '/home/mshen/research/' + sys.argv[1] + '/'
  analyze_contigs(input_folder)


def analyze_contigs(fold):
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


  lengths = []
  num_jumps = 0
  num_single_contigs = 0
  total_num = 0
  single_lengths = []
  covered = [0] * 4700000
  for fil in os.listdir(fold):
    if fnmatch.fnmatch(fil, '*results*'):
      total_num += 1
      with open(fold + fil) as f:
        lines = f.readlines()

      beg = int(lines[0].split()[6])
      end = int(lines[0].split()[7])
      contigs = [[beg, end, 1]]

      for lin in lines[1:]:
        beg = int(lin.split()[6])
        end = int(lin.split()[7])
        if within(contigs[-1][0], contigs[-1][1], beg, end):
          contigs[-1] = expand(contigs[-1][0], contigs[-1][1], beg, end) + [contigs[-1][2] + 1]
        else:
          contigs.append([beg, end, 1])

      print fil, '\n  ', contigs

      best = sorted(contigs, key = lambda l: l[2], reverse = True)[0]
      print '  ', best[1] - best[0], best
      lengths.append(best[1] - best[0])
      num_jumps += len(contigs) - 1
      if len(contigs) == 1:
        num_single_contigs += 1
        single_lengths.append(best[1] - best[0])
      for i in range(best[0], best[1]):
        covered[i] = 1

  print '\nAvg. contig len:', numpy.mean(lengths), ' Stdev:', numpy.std(lengths)
  print 'Single contig avg len:', numpy.mean(single_lengths), ' Stdev:', numpy.std(single_lengths)
  print 'Num. jumps:', num_jumps
  print 'Num. single contigs:', num_single_contigs
  print 'Total num. contigs:', total_num
  print 'Base pairs covered:', sum(covered)



  coverage_outfile = 'out_contig_genome_coverage.out'
  # with open(coverage_outfile, 'w') as f:
  #   f.write('\n'.join([str(s) for s in covered]))

  len_outfile = 'out_contig_lengths_all.out'
  with open(len_outfile, 'w') as f:
    f.write('\n'.join([str(s) for s in lengths]))

  len_outfile_single = 'out_contig_lengths_single.out'
  with open(len_outfile_single, 'w') as f:
    f.write('\n'.join([str(s) for s in single_lengths]))


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start