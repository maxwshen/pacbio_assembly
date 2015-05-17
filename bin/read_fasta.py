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
        curr_read += line.strip()
  if curr_read != '':
    reads.append(curr_read)
  if len(headers) == 0 or len(reads) == 0:
    print 'ERROR: Empty fasta file', fasta_file
  return headers, reads
