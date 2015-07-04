import sys

def read_fasta(fasta_file):
  # Builds header and reads lists from a fasta file
  headers = []
  reads = []
  get_read = False
  with open(fasta_file) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] == '>' or line[0] == '@':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        headers.append(line.strip())
        get_read = True
      elif get_read:
        curr_read += line.strip().upper()
  if curr_read != '':
    reads.append(curr_read)
  if len(headers) == 0 or len(reads) == 0:
    print 'ERROR: Empty fasta file', fasta_file
  return headers, reads

def read_fastq(fastq_fn):
  # There is no official single fastq format, so this may not always work
  headers = []
  reads = []
  quals = []
  get_read = False
  get_qual = False
  with open(fastq_fn) as f:
    curr_read = ''
    curr_qual = ''
    for i, line in enumerate(f):
      if line[0] == '@':
        if curr_qual != '':
          quals.append(curr_qual)
          curr_qual = ''
        headers.append(line.strip())
        get_qual = False
        get_read = True
      elif line == '+':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        get_read = False
        get_qual = True
      else:
        if get_read:
          curr_read += line.strip().upper()
        if get_qual:
          curr_qual += line.strip()
  if curr_read != '':
    reads.append(curr_read)
  if curr_qual != '':
    quals.append(curr_qual)
  if len(headers) == 0 or len(reads) == 0 or len(quals) == 0:
    print 'ERROR: Empty fasta file', fasta_file
  return headers, reads, quals

def reverse_complement(inp):
  if type(inp) != str:
    print 'mylib.reverse_complement Input not a string'
    return ''
  cor = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
  return ''.join([cor[c] for c in inp[::-1]])