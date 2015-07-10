# Perform nhood stats on a header
import itec4
import read_fasta as rf
import datetime

reads_file = '/home/mshen/research/data/reads.20k.rc.fasta'
true_nhoods = '/home/mshen/research/data/genome_nhood_7k.txt'
ktmer_headers_file = '/home/mshen/research/data/20k_v2/temp_ktmer_headersrx_27_6_rc_v2.out'
creads_file = '/home/mshen/research/data/20k_v2/temp_creads.outrx_27_6_rc_v2.out'
# ktmer_headers_file = '/home/yu/max/research/data/20k_v2/temp_ktmer_headersrx_25_4_rc_v2.out'
# creads_file = '/home/yu/max/research/data/20k_v2/temp_creads.outrx_25_4_rc_v2.out'

# ktmer_headers_file = '/home/yu/max/research/data/20k_v2/temp_ktmer_headers55x_15_6_rc_v2_1.out'
# creads_file = '/home/yu/max/research/data/20k_v2/temp_creads.out55x_15_6_rc_v2_1.out'

hr, rr = rf.read_fasta(reads_file)
print datetime.datetime.now()
headers = itec4.build_headers_dict(ktmer_headers_file)
print datetime.datetime.now()
creads = itec4.build_creads_dict(creads_file, hr, rr)
print datetime.datetime.now()
for i in range(len(hr)):
  hr[i] = hr[i].split()[0]
ktmers = headers.keys()

# header = '>R_m140909_094518_42139_c100697390480000001823143803261574_s1_p0/150325/0_32327'

def nhood_stats(base_index, nhood_indices):
  if len(nhood_indices) == 0:
    print 'empty nhood'
    return None
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
    return None
  intersect = len(set([str(s) for s in nhood_indices]).intersection(rest))
  false_reads = [s for s in nhood_indices if str(s) not in rest]
  sensitivity.append(float(intersect) / float(len(rest)))
  specificity.append(float(intersect) / float(len(nhood_indices)))

  numfalse = len(false_reads)
  numMiss = len(rest) - intersect
  print 'false:', len(false_reads)
  print 'missing:', len(rest) - intersect
  return numfalse, numMiss

totals = []
falses = []
misses = []

for i in range(len(hr)):
  if i % 2 == 0:
    header = hr[i]
    nhood, windows = itec4.get_special_1_deg_nhood(header, creads, headers, hr)
    nhood = list(nhood)
    nhood_indices = [hr.index(s) for s in nhood]

    ans = nhood_stats(hr.index(header), nhood_indices)
    if ans is not None:
        (numf, numM) = ans
        totals.append(len(nhood))
        falses.append(numf)
        misses.append(numM)

        print float(sum(totals)) / float(len(totals)), float(sum(falses)) / float(len(falses)), float(sum(misses)) / float(len(misses))

    # nhood = itec4.get_nhood(header, headers, creads, hr)
    # nhood = list(nhood)
    # nhood_indices = [hr.index(s) for s in nhood]
    # itec4.nhood_stats(hr.index(header), nhood_indices)
    # print len(nhood)
    print datetime.datetime.now()
