import sys, string, datetime, random, copy, os, fnmatch, commands
import numpy as np
import read_fasta as rf

def main():
  fold = sys.argv[1]
  files = [s for s in os.listdir(fold) if fnmatch.fnmatch(s, '*hood.fasta')]

  batch_ec(fold)

  # for i in range(len(files)):
  #   print i, len(files)
  #   fn = files[i]
  #   reducer(fold, fn)

def batch_ec(fold):
  # ex reduced nhood file: 40x_999_hood_reduced.fasta 
  # ex base file: 4101_base.fasta
  # NEEDS TO BE RUN INSIDE 1DEG_NHOOD FOLDER
  # BECAUSE yu's EC_tool can't handle directory input
  ec_tool = '/home/yu/program/error_correction_0421.sh'
  blasr_exe = 'blasr'
  blasr_zero = 4      # 0 on debruijn, 4 on Yu's computer
  blasr_zero_len = 8  # 0 on debruijn, 8 on Yu's computer
  blasr_options = '-bestn 1 -m 1 -maxMatch 20'   # Concise output
  leniency = 15     # If more than 15 bp didn't align to genome, error

  reduced_files = [s for s in os.listdir(fold) if fnmatch.fnmatch(s, '*reduced.fasta')]
  for rf in reduced_files:
    num = rf.split('_')[1]
    base = num + '_base.fasta'
    print base, rf
    status = commands.getstatusoutput(ec_tool + ' ' + base + ' ' + rf)[1]
    ec_out = 'C0421_' + base
    status = commands.getstatusoutput('ls ' + ec_out)
    if 'No such file or directory' in status[1]:
      print status
      continue
    new_ec_out = rf.split('.')[0] + '_corr_0421.fasta'
    commands.getstatusoutput('mv ' + ec_out + ' ' + new_ec_out)

    # First with canonical
    e_coli_genome = '/home/yu/e_coli_genome.fasta'     # Canonical
    status = commands.getstatusoutput(blasr_exe + ' ' + new_ec_out + ' ' + e_coli_genome + ' ' + blasr_options)[1]

    if len(status.split()) > blasr_zero_len:
      print status
      acc = status.split()[blasr_zero + 5]
      beg_align = int(status.split()[blasr_zero + 9])
      end_align = int(status.split()[blasr_zero + 10])
      total_len = int(status.split()[blasr_zero + 11])
      if beg_align > leniency or total_len - end_align > leniency:
        print 'canonical failed leniency (15bp) test:' + new_ec_out 
      else:
        print 'canonical ' + num + ': ' + acc + ' ' + new_ec_out

    # Then with quiver
    e_coli_genome = '/home/yu/data/ecoli_consensus_mark.fasta'
    status = commands.getstatusoutput(blasr_exe + ' ' + new_ec_out + ' ' + e_coli_genome + ' ' + blasr_options)[1]

    if len(status.split()) > blasr_zero_len:
      print status
      acc = status.split()[blasr_zero + 5]
      beg_align = int(status.split()[blasr_zero + 9])
      end_align = int(status.split()[blasr_zero + 10])
      total_len = int(status.split()[blasr_zero + 11])
      if beg_align > leniency or total_len - end_align > leniency:
        print 'quiver failed leniency (15bp) test:' + new_ec_out 
      else:
        print 'quiver ' + num + ': ' + acc + ' ' + new_ec_out

    

def reducer(fold, fn):
  step = 5

  with open(fold + fn) as f:
    h, r = rf.read_fasta(fold + fn)
  if len(r) == 0:
    return

  combined = zip(h, r)
  random.shuffle(combined)
  h[:], r[:] = zip(*combined)

  print len(r)

  cutoff = -(len(r) % step)
  if cutoff != 0:
    r = r[:-(len(r) % step)]   # Make r's len divisible by step
    h = h[:-(len(h) % step)]
  print len(r)

  while len(r) > 0:
    print len(r)
    new_f = fold + str(len(r)) + 'x_' + fn.split('.')[0] + '_reduced.fasta'
    with open(new_f, 'w') as f:
      for i in range(len(r)):
        f.write(h[i] + '\n')
        f.write(r[i] + '\n')
    r = r[:-step]
    h = h[:-step]
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
