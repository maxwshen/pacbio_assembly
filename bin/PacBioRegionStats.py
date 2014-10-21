# -*- coding: utf-8 -*-
"""
Created on Sat May 03 18:37:57 2014

@author: Jeffrey
"""

import sys, time

def getGenomicRegion(m5_file_str, size, center, min_overlap):
    '''Takes a blasr -m 5 output file and a region defined by
    the center genomic index of the region (forward strand), and the size
    of the region, (stretches size/2 in both directions from the center).
    Outputs a dictionary of start and end indices of reads that have been aligned
    to this region, cut to only the parts that map to this region.'''
    read_dict = {}
    begin_index = center-size//2
    end_index = center+size//2
    m5_file = open(m5_file_str,'r')
    
    #For debugging
    #genome = readGenome('e_coli_genome.fasta')
    #headers,reads = readFASTA('PacBio_10kb_CLR_mapped_whole_removed_homopolymers.fasta')
    
    for i,line in enumerate(m5_file):
        parts = line.strip().split(" ")
        # read_header = parts[0]
        read_header = '/'.join(parts[0].split('/')[:-1])
        read_length = int(parts[1])
        read_start = int(parts[2])
        read_end = int(parts[3])
        read_strand = parts[4]
        genome_length = int(parts[7])
        genome_start = int(parts[8])
        genome_end = int(parts[9])
        genome_strand = parts[10]
        read_align = parts[17]
        genome_align = parts[19]
        align = parts[18]
        
        flipped = False
        #If on the negative strand of the genome, take the reverse
        #complement of both sequences because genomic index must be forward
        #NOTE: read indices correspond to forward direction. The substring of
        #the read should be taken first before taking the revComp of it.
        if genome_strand == '-':
            
            if read_strand == '+':
                read_strand = '-'
            elif read_strand == '-':
                read_strand = '+'
            genome_strand = '+'
            #Reverse alignment strings
            read_align = revComp(read_align)
            genome_align = revComp(genome_align)
            align = align[::-1]
            
            flipped = True
            
        if genome_end-min_overlap > begin_index and genome_start+min_overlap < end_index:
            '''Given an alignment, a begin_index, and start_indices for the
            parts of the alignment, output the coordinates for the beginning and start
            indices of the read to take.'''
            if not flipped:
                g_i = genome_start  #curr genome index
                r_i = read_start    #curr read index
                a_i = 0             #curr alignment index
                while g_i < begin_index:
                    if read_align[a_i] != '-':
                        r_i += 1
                    if genome_align[a_i] != '-':
                        g_i += 1
                    a_i += 1
                out_read_start = r_i
                a_start = a_i
                g_start = g_i
                
                while g_i < end_index and a_i < len(genome_align):
                    #print '!!!!',end_index,g_i,a_i
                    if read_align[a_i] != '-':
                        r_i += 1
                    if genome_align[a_i] != '-':
                        g_i += 1
                    a_i += 1
                out_read_end = r_i
                a_end = a_i
                g_end = g_i
            else:
                g_i = genome_start
                r_i = read_end
                a_i = 0
                while g_i < begin_index:
                    if read_align[a_i] != '-':
                        r_i -= 1
                    if genome_align[a_i] != '-':
                        g_i += 1
                    a_i += 1
                out_read_end = r_i
                a_start = a_i
                g_start = g_i
                
                while g_i < end_index and a_i < len(genome_align):
                    if read_align[a_i] != '-':
                        r_i -= 1
                    if genome_align[a_i] != '-':
                        g_i += 1
                    a_i += 1
                out_read_start = r_i
                a_end = a_i
                g_end = g_i
                

            # # Debugging
            # print 'Genome, %d:%d' % (g_start,g_end)
            # print genome[g_start:g_end]
            # print 'Read %d, %s, %d:%d' % (headers.index('>'+read_header),read_strand,out_read_start,out_read_end)
            # if read_strand == '+':
            #     print reads[headers.index('>'+read_header)][out_read_start:out_read_end]
            # else:
            #     print revComp(reads[headers.index('>'+read_header)][out_read_start:out_read_end])
            # print
            # print 'Alignment'
            # print genome_align[a_start:a_end]
            # print align[a_start:a_end]
            # print read_align[a_start:a_end]
            # print
            
            read_dict[read_header] = (out_read_start,out_read_end,read_strand,flipped,g_start,g_end)
        
    m5_file.close()
    return read_dict

def extractReads(read_dict, fasta_file_str, out_file_str, write_out):
    fasta_file = open(fasta_file_str,'r')
    if write_out:
        out_file = open(out_file_str,'w')
    currI = (0,0,'+',False,0,0)
    found_header = False
    headers = []
    reads = []
    for i,line in enumerate(fasta_file):
        line = line.strip()
        if line[0] == ">":
            # header = line[1:]
            header = '/'.join(line[1:].split('/')[:-1])
            if header in read_dict:
                found_header = True
                currI = read_dict[header]
            else:
                found_header = False
        elif found_header:
            #print header,currI
            h = ">%s/ex_%d_%d/gi_$%d$_$%d$/" % (header,currI[0],currI[1],currI[4],currI[5])
            if currI[2] == '+':
                r = line[currI[0]:currI[1]]
            elif currI[2] == '-':
                r = revComp(line[currI[0]:currI[1]])

            headers.append(h)
            reads.append(r)
            if write_out:
                out_file.write("%s\n%s\n" % (h,r))
    fasta_file.close()
    if write_out:
        out_file.close()
    return (headers,reads)

def extractGenome(genome_file_str,size,center,out_file_str,write_out,header):
    begin_index = center-size//2
    end_index = center+size//2
    curr_index = 0
    genome_file = open(genome_file_str,'r')
    h = genome_file.readline().strip()
    if not header:
        header = h
    genome = ""
    line = genome_file.readline().strip()
    while line and line != "":
        if curr_index < begin_index and curr_index + len(line) >= begin_index:
            start = begin_index - curr_index
            end = len(line)
            if curr_index + len(line) > end_index:
                end = end_index - curr_index
            genome += line[start:end]
        elif curr_index >= begin_index and curr_index + len(line) <= end_index:
            genome += line
        elif curr_index <= end_index and curr_index + len(line) > end_index:
            end = end_index - curr_index
            genome += line[:end]
        curr_index += len(line)
        line = genome_file.readline().strip()
    genome_file.close()
    if write_out:
        out_file = open(out_file_str,'w')
        out_file.write('%s\n%s\n' % (header,genome))
        out_file.close()
    return genome

def revComp(seq):
    '''Given a DNA sequence, returns the reverse complement of it'''
    trans = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
    out = ''
    for s in seq[::-1]:
        out += trans[s]
    return out

def readGenome(file_name):
    '''Given a fasta file with a genome sequence, returns the genome sequence'''
    g_file = open(file_name,'r')
    header = g_file.readline().strip()
    print header
    genome = ""
    line = g_file.readline().strip()
    while line and line != "":
        genome += line
        line = g_file.readline().strip()
    g_file.close()
    return genome

def readFASTA(file_name):
    '''Given a FASTA file of reads where each sequence is ON ONE LINE,
    return a list of the headers and a list of the reads of the sequences'''
    seq_file = open(file_name,'r')
    headers = []
    seqs = []
    line = seq_file.readline().strip()
    while line and line != "":
        if line[0] == ">":
            headers.append(line)
        else:
            seqs.append(line)
        line = seq_file.readline().strip()
    seq_file.close()
    return (headers,seqs)

def readCenters(file_name):
    c_file = open(file_name,'r')
    centers = []
    kmers = []
    for i,line in enumerate(c_file):
        parts = line.split('_')
        centers.append(int(parts[1][2:]))
        kmers.append(parts[2].split('.')[0])
    c_file.close()
    return centers,kmers

def writeReads(headers,reads,file_name):
    '''Given a list of headers, reads, and a file to write to, writes
    the reads to the file in FASTA format'''
    with open(file_name,'w') as out_file:        
        for i in range(len(headers)):
            gt = ''
            if headers[i][0] != '>':
                gt = '>'
            out_file.write("%s%s\n%s\n" % (gt,headers[i],reads[i]))

def testKmer(genome,kmer):
    for i in range(len(genome)):
        if genome[i:i+len(kmer)] == kmer:
            print i
    
#############################################################################################
    
def main():
    center = int(sys.argv[1])
    size = int(sys.argv[2])
    # min_overlap = size - 100
    min_overlap = 1000
    read_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_removed_homopolymers.fasta'
    genome_file = '/home/mshen/research/data/e_coli_genome.fasta'

    extract(center, size, min_overlap, read_file, genome_file)
    return

def extract(center, size, min_overlap, read_file, genome_file):
    m5_file = '/home/jeyuan/blasr_orig_all_m5.out'
    read_out = 'extracted_reads_c%d_s%d.fasta' % (center,size)
    genome_out = 'extracted_genome_c%d_s%d.fasta' % (center,size)
    write_out = True
    
    fold = '/home/mshen/research/extracts_100k/'
    header = '>ec_genome_region_c%d_s%d/ex_%d_%d/' % (center,size,center-(size//2),center+(size//2))

    #Extract the genome and reads mapping to the genome at a center and size
    read_dict = getGenomicRegion(m5_file, size, center, min_overlap)
    # print read_dict
    headers,reads = extractReads(read_dict, read_file, fold + read_out, write_out)
    genome = extractGenome(genome_file, size, center, fold + genome_out, write_out, header)
    
    return fold + read_out, fold + genome_out

if __name__ == "__main__":
    start_time = time.time()
    main()
    print 'Total time',time.time()-start_time