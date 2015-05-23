import read_fasta as rf

headers_file = '/home/mshen/sample_incorrect_nhoods.txt'
reads_file = '/home/mshen/research/data/reads.20k.rc.fasta'
hr, rr = rf.read_fasta(reads_file)
for i in range(len(hr)):
  hr[i] = hr[i].split()[0]

hs = []
with open(headers_file) as f:
  for i, line in enumerate(f):
    hs.append(line.split())

for i in range(len(hs)):
  h = hs[i]
  out_file = '/home/mshen/sample_incorrect_nhoods_' + str(i) + '.fasta'
  with open(out_file, 'w') as f:
    for h2 in h:
      f.write(h2 + '\n' + rr[hr.index(h2)] + '\n')
