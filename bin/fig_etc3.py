# Finds data on the number of kt-mers per read

import sys
import itec4
import read_fasta as rf

reads_file = '/home/mshen/research/data/reads.20k.rc.fasta'
creads_file = '/home/mshen/research/data/20k_v2/temp_creads.outrx_27_6_rc_v2.out'
hr, rr = rf.read_fasta(reads_file)
creads_dict = itec4.build_creads_dict(creads_file, hr, rr)
for i in range(len(hr)):
  hr[i] = hr[i].split()[0]
for i in range(20):
  print hr[i], len(creads_dict[hr[i]]) - 1 / 2, len(rr[i])

sys.exit(0)

data = []
for i in range(len(creads_dict.keys())):
  k = creads_dict.keys()[i]
  data.append((len(creads_dict[k]) - 1) / 2)
  if i % 5000 == 0:
    print i

print '\n'.join([str(s) for s in data])
print float(sum(data)) / float(len(data))