# Perform nhood stats on a header
import itec4
import read_fasta as rf
import datetime

reads_file = '/home/yu/max/research/data/reads.20k.rc.fasta'
true_nhoods = '/home/yu/max/research/data/genome_nhood.txt'
# ktmer_headers_file = '/home/yu/max/research/data/20k_v2/temp_ktmer_headersrx_27_6_rc_v2.out'
# creads_file = '/home/yu/max/research/data/20k_v2/temp_creads.outrx_27_6_rc_v2.out'
ktmer_headers_file = '/home/yu/max/research/data/20k_v2/temp_ktmer_headersrx_25_4_rc_v2.out'
creads_file = '/home/yu/max/research/data/20k_v2/temp_creads.outrx_25_4_rc_v2.out'

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

for i in range(len(hr)):
  if i % 2 == 0:
    header = hr[i]
    nhood, windows = itec4.get_simple_1_deg_nhood(header, creads, headers, hr)
    nhood = list(nhood)
    nhood_indices = [hr.index(s) for s in nhood]

    itec4.nhood_stats(hr.index(header), nhood_indices)
    print len(nhood)

    nhood = itec4.get_nhood(header, headers, creads, hr)
    nhood = list(nhood)
    nhood_indices = [hr.index(s) for s in nhood]
    itec4.nhood_stats(hr.index(header), nhood_indices)
    print len(nhood)
    print datetime.datetime.now()
