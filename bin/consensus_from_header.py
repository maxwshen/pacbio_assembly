import sys, string, datetime, random, copy, os, commands
import itec4
import read_fasta as rf

def main():
  header = '>' + sys.argv[1]
  e_coli_genome = '/home/mshen/research/data/e_coli_genome.fasta'
  # ec_tool = '/home/mshen/research/bin/error_correction_1218.sh'
  # reads_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
  # creads_file = '/home/mshen/research/data/22.4_creads.out'
  # ktmer_headers_file = '/home/mshen/research/data/22.4_ktmer_headers.out'
  blasr_exe = '/home/jeyuan/blasr/alignment/bin/blasr'
  blasr_options = '-bestn 1 -m 1'   # Concise output
  temp_sig = str(datetime.datetime.now()).split()[1]

  # New dataset
  ec_tool = '/home/lin/program/error_correction_5X_0204.sh'
  reads_file = '/home/mchaisso/datasets/pacbio_ecoli/reads.20k.fasta'
  # creads_file = '/home/mshen/research/data/22.8_creads_20k.out'
  # ktmer_headers_file = '/home/mshen/research/data/22.8_ktmer_headers_20k.out'
  creads_file = '/home/mshen/research/data/temp_creads.out_28_6.out'
  ktmer_headers_file = '/home/mshen/research/data/temp_ktmer_headers_28_6.out'

  creads = itec4.build_creads_dict(creads_file, reads_file)
  headers = itec4.build_headers_dict(ktmer_headers_file)
  hr, rr = rf.read_fasta(reads_file)
  # Compensate for new dataset
  for i in range(len(hr)):
    hr[i] = hr[i].split()[0]

  con = itec4.error_correct(ec_tool, header, headers, creads, hr, rr, temp_sig_out = temp_sig)
  if len(con) == 0:
    print 'FAILURE IN ERROR CORRECTION'
    sys.exit(0)

  return

  temp_file = 'temp_cfh_' + temp_sig + '.fasta'
  temp2_file = 'temp_cfh2_' + temp_sig + '.fasta'
  with open(temp_file, 'w') as f:
    f.write(header + '\n' + con)

  status = commands.getstatusoutput(blasr_exe + ' ' + temp_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]
  if len(status) != 0:
    print status

  collected_h = set()
  ktmers = []
  if header not in creads or len(creads[header]) == 1:
    pass
  for i in range(len(creads[header])):
    if i % 2 == 1:
      ktmers.append(creads[header][i])
  for kt in ktmers:
    for h in headers[kt]:
      collected_h.add(h)

  to_con = []
  to_gen = []
  for ch in collected_h:
    with open(temp2_file, 'w') as f:
      f.write(ch + '\n' + rr[hr.index(ch)])
    status = commands.getstatusoutput(blasr_exe + ' ' + temp2_file + ' ' + temp_file + ' ' + blasr_options)[1]
    to_con.append(status)
    status = commands.getstatusoutput(blasr_exe + ' ' + temp2_file + ' ' + e_coli_genome + ' ' + blasr_options)[1]
    to_gen.append(status)

  print sum([1 for s in to_con if len(s) > 0]), 'used in consensus out of', len(to_con)
  for tg in to_gen:
    print tg


if __name__ == '__main__':
  # Initiates program and records total time
  start = datetime.datetime.now()
  # print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  # print '\n\nEnd:', end, '\nTotal:', end - start