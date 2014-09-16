# -*- coding: utf-8 -*-
"""
Created on Sat May 03 18:37:57 2014

@author: Jeffrey
"""

import time,gc

def getGenomicRegion(m5_file_str,size,center,min_overlap):
    '''Takes a blasr -m 5 output file and a region defined by
    the center genomic index of the region - forward strand, and the size
    of the region, (stretches size/2 in both directions from the center).
    Outputs a dictionary of start and end indices of reads that have been aligned
    to this region, cut to only the parts that map to this region.'''
    reads = {}
    begin_index = center-size//2
    end_index = center+size//2
    m5_file = open(m5_file_str,'r')
    for i,line in enumerate(m5_file):
        parts = line.strip().split(" ")
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
        #align = parts[18]
        
        flipped = False
        #If on the negative strand of the genome, take the reverse
        #complement of both sequences
        if genome_strand == '-':
            '''gr_start = genome_length - genome_end
            gr_end = genome_length - genome_start
            genome_start = gr_start
            genome_end = gr_end'''
            rr_start = read_length - read_end
            rr_end = read_length - read_start
            read_start = rr_start
            read_end = rr_end
            
            if read_strand == '+':
                read_strand = '-'
            elif read_strand == '-':
                read_strand = '+'
            #Reverse alignment strings
            read_align = revComp(read_align)
            genome_align = revComp(genome_align)
            #align = align[::-1]
            
            flipped = True
            
        if genome_end-min_overlap > begin_index and genome_start+min_overlap < end_index:
            '''Given an alignment, a begin_index, and start_indices for the
            parts of the alignment, output the coordinates for the beginning and start
            indices of the read to take.'''
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
            #a_start = a_i
            
            while g_i < end_index and a_i < len(genome_align):
                #print '!!!!',end_index,g_i,a_i
                if read_align[a_i] != '-':
                    r_i += 1
                if genome_align[a_i] != '-':
                    g_i += 1
                a_i += 1
            out_read_end = r_i
            #a_end = a_i

            
            '''#Debugging
            print out_read_start,read_start
            print out_read_end,read_end
            print genome_strand,read_strand
            print read_align[a_start:a_end]
            print align[a_start:a_end]
            print genome_align[a_start:a_end]
            print len(read_align[a_start:a_end])
            print len(''.join(read_align[a_start:a_end].split('-')))
            print len(genome_align[a_start:a_end])
            print len(''.join(genome_align[a_start:a_end].split('-')))
            print'''
            reads[read_header] = (out_read_start,out_read_end,read_strand,flipped)
        
    m5_file.close()
    return reads

def extractReads(read_dict,fasta_file_str,out_file_str,write_out):
    fasta_file = open(fasta_file_str,'r')
    if write_out:
        out_file = open(out_file_str,'w')
    currI = (0,0,'+',True)
    found_header = False
    headers = []
    reads = []
    for i,line in enumerate(fasta_file):
        line = line.strip()
        if line[0] == ">":
            header = line[1:]
            if header in read_dict:
                found_header = True
                currI = read_dict[header]
            else:
                found_header = False
        elif found_header:
            #print header,currI
            h = ">%s/%d_%d/" % (header,currI[0],currI[1])
            if currI[2] == '+':
                r = line[currI[0]:currI[1]]
            elif currI[2] == '-':
                r = revComp(line)[currI[0]:currI[1]]

            headers.append(h)
            reads.append(r)
            if write_out:
                out_file.write("%s\n%s\n" % (h,r))
            #print h
            #print r
    fasta_file.close()
    if write_out:
        out_file.close()
    return (headers,reads)

def extractGenome(genome_file_str,size,center,out_file_str,write_out,header):
    begin_index = center-size//2
    end_index = center+size//2
    curr_index = 0
    genome_file = open(genome_file_str,'r')
    #header = genome_file.readline().strip()
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

def findKTmers(k,t,reads,genome,out_file_str):
    '''Given a k, t, a list of headers and reads, and a genome, find
    the genomic index of the starting position of all the (k,t)-mers
    that exactly match to the genome'''
    allkmers = {}
    for j in range(len(reads)):
        for i in range(len(reads[j])-k+1):
            kmer = reads[j][i:i+k]
            if not allkmers.has_key(kmer):
                allkmers[kmer] = [1,0]
            else:
                allkmers[kmer][0] += 1
        #if j % 1000 == 0:
        #    print '%d out of %d reads' % (j,len(reads))
    print 'Done reading reads'
    out_file = open(out_file_str,'w')
    out_file.write("Genomic positions of (%d,%d)-mers:\n\n" % (k,t))
    for i in range(len(genome)-k+1):
        kmer = genome[i:i+k]
        if kmer in allkmers and allkmers[kmer][0] >= t:
            allkmers[kmer][1] += 1
            out_file.write('%s, t=%d, n=%d time in genome: %d\n' % \
                        (kmer,allkmers[kmer][0],allkmers[kmer][1],i))
            out_file.flush()
        #if i % 1000 == 0:
        #    print 'genome',i
    out_file.close()
                
def checkHomopolymers(genome):
    '''Given a genome, count and print the longest stretch of homopolymer
    repeats for each of the 4 bases'''
    currLength = 0
    lens = {'A':0,'C':0,'G':0,'T':0}
    currBase = ''
    for c in genome:
        if c == currBase:
            currLength += 1
        else:
            if currBase:
                if currLength > lens[currBase]:
                    lens[currBase] = currLength
                currLength = 0
            currBase = c
    if currLength > lens[currBase]:
        lens[currBase] = currLength
    for c in lens:
        print c,lens[c]

def genDelKmers(kmer,del_size):
    '''Generate a list of all kmers made from del <= del_size, where
    del_kmers[0] = no deletions in k
    del_kmers[1] = 1 deletion in k ... up to del_kmers[del_size]'''
    k = len(kmer)
    del_kmers = [[] for x in range(del_size+1)]
    del_kmers[0] = [kmer]
    next_kmers = []
    for j in range(del_size):
        for curr in del_kmers[j]:
            for i in range(1,k-j-1):
                next_kmers.append(curr[:i]+curr[i+1:])
        del_kmers[j+1] = set(next_kmers[:])
        next_kmers = []
    return del_kmers

def genInsKmers(kmer,ins_size):
    '''Generate a list of all kmers made from ins <= ins_size, where
    ins_kmers[0] = no insertions in k
    ins_kmers[1] = 1 insertion in k ... up to ins_kmers[ins_size]'''
    k = len(kmer)
    alpha = ['A','C','G','T']
    ins_kmers = [[] for x in range(ins_size+1)]
    ins_kmers[0] = [kmer]
    next_kmers = []
    for j in range(ins_size):
        for curr in ins_kmers[j]:
            for i in range(1,k-j-1):
                for m in range(len(alpha)):
                    next_kmers.append(curr[:i]+alpha[m]+curr[i:])
        ins_kmers[j+1] = set(next_kmers[:])
        next_kmers = []
    return ins_kmers

def genSubKmers(kmer,new_k):
    '''Generate a list of all "sub kmers" made from the original
    kmer of size new_k'''
    sub_kmers = set([kmer[x:x+new_k] for x in range(len(kmer)-new_k+1)])
    return sub_kmers

def delCheck(reads,genome,k,del_size,offset,out_file_str):
    '''Takes a set of reads mapped to a 1000 bp region of the genome
    centered at a kmer at the 500 bp position, a k, the number of
    deletions to kmers to consider, an offset for
    the tolerance for mapping to positions, and an ouput filename.
    This function finds all of the kmers formed by k = 0, 1, ... del_size
    deletions to the central k-mer, checks whether they show up in the reads,
    and prints and writes summary statistics for whether the deletion k-mers
    in the reads map to the correct location of the genome'''
    out_file = open(out_file_str,'a')
    g_mid = int(len(genome)/2)
    kmer = genome[g_mid:g_mid+k]
    out_file.write('%s\n' % kmer)
    print kmer
    
    curr_kmers = genDelKmers(kmer,del_size)
    
    success = [0 for x in range(del_size+1)]
    failure = [0 for x in range(del_size+1)]
    for i in range(len(reads)):
        for j in range(len(reads[i])-k+1+del_size):
            found_short = False
            for m in range(del_size+1):
                if j < len(reads[i])-k+1+m and reads[i][j:j+k-m] in curr_kmers[m]:
                    out_file.write('%s del=%d, %d of %d ' % (reads[i][j:j+k-m],m,j,len(reads[i])))
                    print reads[i][j:j+k-m],'del=%d,'%m,j,'of',len(reads[i]),
                    if abs(g_mid - j) <= abs(len(reads[i])-len(genome))+offset:
                        success[m] += 1
                        out_file.write('S \n')
                        print 'S'
                    else:
                        failure[m] += 1   
                        out_file.write('F \n')
                        print 'F'
                    found_short = True
                    break
            if found_short:
                break

                    
    out_file.write('del =\t\t')
    print 'del =\t',
    for i in range(len(success)):
        out_file.write('%d\t' % i)
        print '%d\t' % i,
    
    out_file.write('\nSuccesses\t')
    print '\nSuccesses\t',
    for i in range(len(success)):
        out_file.write('%d\t' % success[i])
        print '%d\t' % success[i],

    out_file.write('\nFailures\t')
    print '\nFailures\t',
    for i in range(len(success)):
        out_file.write('%d\t' % failure[i])
        print '%d\t' % failure[i],
    out_file.write('\n\n')
    print '\n'
    out_file.write('Total successes:\t%d\n' % sum(success))
    out_file.write('Total failures:\t%d\n\n' % sum(failure))
    print 'Total num successes:\t%d' % sum(success)
    print 'Total num failures:\t%d\n' % sum(failure)
    out_file.close()
    
def indelCheck(reads,genome,k,indel_size,offset,out_file_str,delOrIns):
    '''Takes a set of reads mapped to a 1000 bp region of the genome
    centered at a kmer at the 500 bp position, a k, the number of
    insertions to kmers to consider, an offset for
    the tolerance for mapping to positions, and an ouput filename.
    This function finds all of the kmers formed by k = 0, 1, ... ins_size
    insertions to the central k-mer, checks whether they show up in the reads,
    and prints and writes summary statistics for whether the insertion k-mers
    in the reads map to the correct location of the genome
    If delOrIns is true, then do delCheck, if false then do insCheck'''
    out_file = open(out_file_str,'a')
    g_mid = int(len(genome)/2)
    kmer = genome[g_mid:g_mid+k]
    out_file.write('%s\n' % kmer)
    print kmer
    
    if delOrIns:
        curr_kmers = genDelKmers(kmer,indel_size)
        adj_k = k-indel_size
        out_str = 'del'
    else:
        curr_kmers = genInsKmers(kmer,indel_size)
        adj_k = k+indel_size
        out_str = 'ins'
    
    success = [0 for x in range(indel_size+1)]
    failure = [0 for x in range(indel_size+1)]
    for i in range(len(reads)):
        for j in range(len(reads[i])-adj_k+1):
            found_short = False
            for m in range(indel_size+1):
                if delOrIns:
                    adj_m = m
                else:
                    adj_m = -1*m
                if j < len(reads[i])-k+1+adj_m and reads[i][j:j+k-adj_m] in curr_kmers[m]:
                    out_file.write('%s %s=%d, %d of %d ' % (reads[i][j:j+k-adj_m],out_str,m,j,len(reads[i])))
                    print reads[i][j:j+k-adj_m],'%s=%d,'%(out_str,m),j,'of',len(reads[i]),
                    if abs(g_mid - j) <= abs(len(reads[i])-len(genome))+offset:
                        success[m] += 1
                        out_file.write('S \n')
                        print 'S'
                    else:
                        failure[m] += 1   
                        out_file.write('F \n')
                        print 'F'
                    found_short = True
                    break
            if found_short:
                break

                    
    out_file.write('%s =\t\t'%out_str)
    print '%s =\t'%out_str,
    for i in range(len(success)):
        out_file.write('%d\t' % i)
        print '%d\t' % i,
    
    out_file.write('\nSuccesses\t')
    print '\nSuccesses\t',
    for i in range(len(success)):
        out_file.write('%d\t' % success[i])
        print '%d\t' % success[i],

    out_file.write('\nFailures\t')
    print '\nFailures\t',
    for i in range(len(success)):
        out_file.write('%d\t' % failure[i])
        print '%d\t' % failure[i],
    out_file.write('\n\n')
    print '\n'
    out_file.write('Total successes:\t%d\n' % sum(success))
    out_file.write('Total failures:\t%d\n\n' % sum(failure))
    print 'Total num successes:\t%d' % sum(success)
    print 'Total num failures:\t%d\n' % sum(failure)
    out_file.close()

def subCheck(reads,genome,k,new_k,offset,out_file_str):
    '''Takes a set of reads mapped to a 1000 bp region of the genome
    centered at a kmer at the 500 bp position, a k, the number of
    insertions to kmers to consider, an offset for
    the tolerance for mapping to positions, and an ouput filename.
    This function finds all of the kmers formed by k = 0, 1, ... ins_size
    insertions to the central k-mer, checks whether they show up in the reads,
    and prints and writes summary statistics for whether the insertion k-mers
    in the reads map to the correct location of the genome'''
    out_file = open(out_file_str,'a')
    g_mid = int(len(genome)/2)
    kmer = genome[g_mid:g_mid+k]
    out_file.write('%s\n' % kmer)
    print kmer
    
    curr_kmers = genSubKmers(kmer,new_k)
    
    success = 0
    failure = 0
    for i in range(len(reads)):
        for j in range(len(reads[i])+new_k-1):
            if j < len(reads[i])+new_k-1 and reads[i][j:j+new_k] in curr_kmers:
                out_file.write('%s, %d of %d ' % (reads[i][j:j+new_k],j,len(reads[i])))
                print reads[i][j:j+new_k],j,'of',len(reads[i]),
                if abs(g_mid - j) <= abs(len(reads[i])-len(genome))+offset:
                    success += 1
                    out_file.write('S \n')
                    print 'S'
                else:
                    failure += 1   
                    out_file.write('F \n')
                    print 'F'
                break

                    
    out_file.write('new_k = %d\nSuccesses = %d\nFailures = %d\n' % (new_k,success,failure))
    print 'new_k = %d\nSuccesses = %d\nFailures = %d' % (new_k,success,failure)
    
    out_file.close()

def delCheckRange(m5_file,read_file,genome_file,gi_file,r_size,r_center,size,min_overlap,k,del_size,offset,out_file_str):
    '''This function enables the calculation of all solid k-mers whether r_size of
    the center k-mer, generates a genome and reads centered on each of those solid k-mers,
    and runs delCheck,insCheck,or subCheck on those centered genome and reads, determining how many
    k-mers with deletions occur in the reads around each solid k-mer.'''
    
    solids = findSolids(gi_file,r_size,r_center)
    str_out = ''
    write_out = False
    with open(out_file_str,'w') as out_file:
        out_file.write('k = %d, num_del = %d\n\n' % (k,del_size))
        print 'k = %d, num_del = %d\n' % (k,del_size)
    for tup in solids:
        with open(out_file_str,'a') as out_file:
            out_file.write("kmer=%s, t=%d, n=%d, gi=%d\n" % tup)
        print "kmer=%s, t=%d, n=%d, gi=%d" % tup
        center = tup[-1]
        read_dict = getGenomicRegion(m5_file,size,center,min_overlap)
        headers,reads = extractReads(read_dict,read_file,str_out,write_out)
        genome = extractGenome(genome_file,size,center,str_out,write_out)
        delCheck(reads,genome,k,del_size,offset,out_file_str)
        #kmer = tup[0]
        #delFix(reads,genome,k,del_size,offset,kmer,500,out_file_str)

def findSolids(gi_file_str,size,center):
    '''Takes a genomic indices file, a center, and a region size to consider
    and outputs the indices of all of the solid ktmers found in this region.'''
    begin_index = center-size//2
    end_index = center+size//2
    gi_file = open(gi_file_str,'r')
    solids = []
    for i,line in enumerate(gi_file):
        if i > 1:
            parts = line.strip().split(' ')
            g_i = int(parts[-1])
            if g_i >= begin_index and g_i < end_index:
                kmer = parts[0][:-1]
                t = int(parts[1][:-1].split('=')[-1])
                n = int(parts[2].split('=')[-1])
                solids.append((kmer,t,n,g_i))
            
    return solids

def resolveFixes(kmer_fix,reads,solids,k):
    '''Given a dictionary of what reads to fix, with the indices, number of
    deletions, and kmer to insert, return a new set of reads with the fixes.'''
    pass
    
def delFixIterate(m5_file,read_file,genome_file,size,center,min_overlap,k,t,del_size,offset,num_iter,out_file_str,out_read_file):
    '''Given an m5_file, a read and a genome file, first extract the region of
    the reads and genome using size, center, min_overlap. Then find all solid
    (k,t)-mers in that region and where they occur, and detect and fix
    occurrences of (k,t)-mers that have up to del_size deletions. Repeat the
    previous step num_iter times.'''
    str_out = ''
    ktmer_out = 'del_ktmer_temp.txt'
    #g_out = 'extracted_genome_c%d_s%d.fasta' % (center,size)
    r_size = size
    r_center = 500
    
    read_dict = getGenomicRegion(m5_file,size,center,min_overlap)
    headers,reads = extractReads(read_dict,read_file,str_out,False)
    genome = extractGenome(genome_file,size,center,str_out,False)
    
    num_solid = [0 for x in range(num_iter)]
    exact_s = [0 for x in range(num_iter)]
    exact_f = [0 for x in range(num_iter)]
    successes = [0 for x in range(num_iter)]
    failures = [0 for x in range(num_iter)]
    
    with open(out_file_str,'w') as out_file:
        pass
    
    for i in range(num_iter):
        findKTmers(k,t,reads,genome,ktmer_out)
        solids = findSolids(ktmer_out,r_size,r_center)
        num_solid[i] = len(solids)
        #kmer_fix = {}

        with open(out_file_str,'a') as out_file:
            out_file.write('k = %d, t = %d, num_del = %d\n' % (k,t,del_size))
            out_file.write('Total ktmers %d, Iteration %d:\n\n' % (num_solid[i],i+1))
            print 'k = %d, t = %d, num_del = %d' % (k,t,del_size)
            print 'Total ktmers %d, Iteration %d:\n' % (num_solid[i],i+1)
        for j in range(len(solids)):
            with open(out_file_str,'a') as out_file:
                out_file.write("kmer=%s, t=%d, n=%d, gi=%d\n" % solids[j])
            print "kmer=%s, t=%d, n=%d, gi=%d" % solids[j]
            kmer = solids[j][0]
            index = solids[j][-1]
            reads,ex_succ,ex_fail,succ,fail = delFix(reads,genome,k,del_size,offset,kmer,index,out_file_str)
            '''for tup in read_changes:
                if tup[0] not in kmer_fix:
                    kmer_fix[tup[0]] = [tup[1:]+[j]]
                else:
                    kmer_fix[tup[0]] += [tup[1:]+[j]]'''
            exact_s[i] += ex_succ
            exact_f[i] += ex_fail
            successes[i] += succ
            failures[i] += fail
        
        with open(out_file_str,'a') as out_file:
            out_file.write('Iter %d:\n' % (i+1))
            out_file.write('Num ktmers %d\n' % num_solid[i])
            out_file.write('Exact %dS + %dF\n' % (exact_s[i],exact_f[i]))
            out_file.write('Success %d\n' % successes[i])
            out_file.write('Failure %d\n\n' % failures[i])
            print 'Iter %d:' % (i+1)
            print 'Num ktmers %d' % num_solid[i]
            print 'Exact %dS + %dF' % (exact_s[i],exact_f[i])
            print 'Success %d' % successes[i]
            print 'Failure %d\n' % failures[i]
        
        #print kmer_fix
        #new_reads = resolveFixes(kmer_fix,reads,solids,k)
        
    with open(out_file_str,'a') as out_file:
        out_file.write('Summary Statistics:\n')
        out_file.write('Iteration\t')
        print 'Iteration\t',
        for i in range(num_iter):
            out_file.write('%d\t' % (i+1))
            print '%d\t' % (i+1),
        out_file.write('Tot\t')
        print 'Tot\t',
            
        out_file.write('\nNum ktmers\t')
        print '\nNum ktmers\t',
        for i in range(num_iter):
            out_file.write('%d\t' % num_solid[i])
            print '%d\t' % num_solid[i],
        out_file.write('%d\t' % max(num_solid))
        print '%d\t' % max(num_solid),
    
        out_file.write('\nExact\t')
        print '\nExact\t',
        for i in range(num_iter):
            out_file.write('%dS+%dF\t' % (exact_s[i],exact_f[i]))
            print '%dS+%dF\t' % (exact_s[i],exact_f[i]),
        out_file.write('%dS+%dF\t' % (max(exact_s),max(exact_f)))
        print '%dS+%dF\t' % (max(exact_s),max(exact_f)),

        out_file.write('\nSuccess\t')
        print '\nSuccess\t',
        for i in range(num_iter):
            out_file.write('%d\t' % successes[i])
            print '%d\t' % successes[i],
        out_file.write('%d\t' % sum(successes))
        print '%d\t' % sum(successes),

        out_file.write('\nFailure\t')
        print '\nFailure\t',
        for i in range(num_iter):
            out_file.write('%d\t' % failures[i])
            print '%d\t' % failures[i],
        out_file.write('%d\t\n\n' % sum(failures))
        print '%d\t\n' % sum(failures)

    writeReads(headers,reads,out_read_file)

def indelFixIterate(m5_file,read_file,genome_file,size,center,min_overlap,k,t,indel_size,offset,num_iter,out_file_str,out_read_file,indel_state,fixIncorrect):
    '''Given an m5_file, a read and a genome file, first extract the region of
    the reads and genome using size, center, min_overlap. Then find all solid
    (k,t)-mers in that region and where they occur, and detect and fix
    occurrences of (k,t)-mers that have up to del_size deletions. Repeat the
    previous step num_iter times. If indel_state is 1, then fix insertions;
    if -1, then fix deletions, if 0, then fix both. If fixIncorrect is true,
    fix k-mers found even at the wrong positions.'''
    str_out = ''
    ktmer_out = 'del_ktmer_temp.txt'
    g_out = 'extracted_genome_c%d_s%d.fasta' % (center,size)
    r_out = 'extracted_reads_c%d_s%d.fasta' % (center,size)
    r_size = size
    r_center = 500
    if indel_state == 1:
        indel_str = "ins"
        delOrIns = False
    elif indel_state == -1:
        indel_str = "del"
        delOrIns = True
    elif indel_state == 0:
        indel_str = 'both'
        delOrIns = True
        
        exact_s_i = [0 for x in range(num_iter)]
        exact_f_i = [0 for x in range(num_iter)]
        successes_i = [0 for x in range(num_iter)]
        failures_i = [0 for x in range(num_iter)]
    else:
        pass
    
    read_dict = getGenomicRegion(m5_file,size,center,min_overlap)
    headers,reads = extractReads(read_dict,read_file,r_out,True)
    genome = extractGenome(genome_file,size,center,g_out,True)
    
    num_solid = [0 for x in range(num_iter)]
    exact_s = [0 for x in range(num_iter)]
    exact_f = [0 for x in range(num_iter)]
    successes = [0 for x in range(num_iter)]
    failures = [0 for x in range(num_iter)]
    
    with open(out_file_str,'w') as out_file:
        pass
    
    for i in range(num_iter):
        findKTmers(k,t,reads,genome,ktmer_out)
        solids = findSolids(ktmer_out,r_size,r_center)
        num_solid[i] = len(solids)

        with open(out_file_str,'a') as out_file:
            out_file.write('k = %d, t = %d, num_indel = %d, state = %s\n' % (k,t,indel_size,indel_str))
            out_file.write('Total ktmers %d, Iteration %d:\n\n' % (num_solid[i],i+1))
            print 'k = %d, t = %d, num_del = %d, state = %s' % (k,t,del_size,indel_str)
            print 'Total ktmers %d, Iteration %d:\n' % (num_solid[i],i+1)
        for j in range(len(solids)):
            with open(out_file_str,'a') as out_file:
                out_file.write("kmer=%s, t=%d, n=%d, gi=%d\n" % solids[j])
            print "kmer=%s, t=%d, n=%d, gi=%d" % solids[j]
            kmer = solids[j][0]
            index = solids[j][-1]
            
            reads,ex_succ,ex_fail,succ,fail = indelFix(reads,genome,k,indel_size,offset,kmer,index,out_file_str,delOrIns,fixIncorrect)
            exact_s[i] += ex_succ
            exact_f[i] += ex_fail
            successes[i] += succ
            failures[i] += fail
            
            if indel_state == 0:
                reads,ex_succ,ex_fail,succ,fail = indelFix(reads,genome,k,indel_size,offset,kmer,index,out_file_str,False,fixIncorrect)
                exact_s_i[i] += ex_succ
                exact_f_i[i] += ex_fail
                successes_i[i] += succ
                failures_i[i] += fail
        
        with open(out_file_str,'a') as out_file:
            out_file.write('Iter %d:\n' % (i+1))
            out_file.write('Num ktmers %d\n' % num_solid[i])
            print 'Iter %d:' % (i+1)
            print 'Num ktmers %d' % num_solid[i]
            
            if indel_state == -1 or indel_state == 0:
                out_file.write('Deletions:\n')
                print 'Deletions:'
            elif indel_state == 1:
                out_file.write('Insertions:\n')
                print 'Insertions:'
                
            out_file.write('Exact successes %d\n' % exact_s[i])
            out_file.write('Exact failures %d\n' % exact_f[i])            
            out_file.write('Successes %d\n' % successes[i])
            out_file.write('Failures %d\n\n' % failures[i])
            
            print 'Exact successes %d' % exact_s[i]
            print 'Exact failures %d' % exact_f[i]
            print 'Successes %d' % successes[i]
            print 'Failures %d\n' % failures[i]
            
            if indel_state == 0:
                out_file.write('Insertions:\n')
                print 'Insertions:'
                
                out_file.write('Exact successes %d\n' % exact_s_i[i])
                out_file.write('Exact failures %d\n' % exact_f_i[i])            
                out_file.write('Success %d\n' % successes_i[i])
                out_file.write('Failure %d\n\n' % failures_i[i])
                
                print 'Exact successes %d' % exact_s_i[i]
                print 'Exact failures %d' % exact_f_i[i]
                print 'Successes %d' % successes_i[i]
                print 'Failures %d\n' % failures_i[i]
                
                out_file.write('Total:\n')
                print 'Total:'
                
                out_file.write('Exact successes %d\n' % (exact_s[i] + exact_s_i[i]))
                out_file.write('Exact failures %d\n' % (exact_f[i] + exact_f_i[i]))
                out_file.write('Success %d\n' % (successes[i] + successes_i[i]))
                out_file.write('Failure %d\n\n' % (failures[i] + failures_i[i]))
                
                print 'Exact successes %d' % (exact_s[i] + exact_s_i[i])
                print 'Exact failures %d' % (exact_f[i] + exact_f_i[i])
                print 'Successes %d' % (successes[i] + successes_i[i])
                print 'Failures %d\n' % (failures[i] + failures_i[i])

        
    with open(out_file_str,'a') as out_file:
        out_file.write('Summary Statistics:\n')
        out_file.write('Iteration\t')
        print 'Iteration\t',
        for i in range(num_iter):
            out_file.write('%d\t' % (i+1))
            print '%d\t' % (i+1),
        out_file.write('Tot\t')
        print 'Tot\t',
            
        out_file.write('\nNum ktmers\t')
        print '\nNum ktmers\t',
        for i in range(num_iter):
            out_file.write('%d\t' % num_solid[i])
            print '%d\t' % num_solid[i],
        out_file.write('%d\t' % max(num_solid))
        print '%d\t' % max(num_solid),
    
        if indel_state == -1 or indel_state == 0:
            out_file.write('\nDeletions:')
            print '\nDeletions:',
        elif indel_state == 1:
            out_file.write('\nInsertions:')
            print '\nInsertions:',
        
        out_file.write('\nExact Succ\t')
        print '\nExact Succ\t',
        for i in range(num_iter):
            out_file.write('%d\t' % exact_s[i])
            print '%d\t' % exact_s[i],
        out_file.write('%d\t' % max(exact_s))
        print '%d\t' % max(exact_s),

        out_file.write('\nExact Fail\t')
        print '\nExact Fail\t',
        for i in range(num_iter):
            out_file.write('%d\t' % exact_f[i])
            print '%d\t' % exact_f[i],
        out_file.write('%d\t' % max(exact_f))
        print '%d\t' % max(exact_f),

        out_file.write('\nSuccess\t')
        print '\nSuccesses  \t',
        for i in range(num_iter):
            out_file.write('%d\t' % successes[i])
            print '%d\t' % successes[i],
        out_file.write('%d\t' % sum(successes))
        print '%d\t' % sum(successes),

        out_file.write('\nFailure\t')
        print '\nFailures  \t',
        for i in range(num_iter):
            out_file.write('%d\t' % failures[i])
            print '%d\t' % failures[i],
        out_file.write('%d\t' % sum(failures))
        print '%d\t' % sum(failures)
        
        if indel_state == 0:
            out_file.write('\nInsertions:')
            print '\nInsertions:',
            
            out_file.write('\nExact Succ\t')
            print '\nExact Succ\t',
            for i in range(num_iter):
                out_file.write('%d\t' % exact_s_i[i])
                print '%d\t' % exact_s_i[i],
            out_file.write('%d\t' % max(exact_s_i))
            print '%d\t' % max(exact_s_i),
    
            out_file.write('\nExact Fail\t')
            print '\nExact Fail\t',
            for i in range(num_iter):
                out_file.write('%d\t' % exact_f_i[i])
                print '%d\t' % exact_f_i[i],
            out_file.write('%d\t' % max(exact_f_i))
            print '%d\t' % max(exact_f_i),
    
            out_file.write('\nSuccess\t')
            print '\nSuccesses  \t',
            for i in range(num_iter):
                out_file.write('%d\t' % successes_i[i])
                print '%d\t' % successes_i[i],
            out_file.write('%d\t' % sum(successes_i))
            print '%d\t' % sum(successes_i),
    
            out_file.write('\nFailure\t')
            print '\nFailures  \t',
            for i in range(num_iter):
                out_file.write('%d\t' % failures_i[i])
                print '%d\t' % failures_i[i],
            out_file.write('%d\t' % sum(failures_i))
            print '%d\t' % sum(failures_i)
            
            
            out_file.write('\nTotal:')
            print '\nTotal:',
            
            out_file.write('\nExact Succ\t')
            print '\nExact Succ\t',
            for i in range(num_iter):
                out_file.write('%d\t' % (exact_s[i]+exact_s_i[i]))
                print '%d\t' % (exact_s[i]+exact_s_i[i]),
            out_file.write('%d\t' % (max(exact_s) + max(exact_s_i)))
            print '%d\t' % (max(exact_s) + max(exact_s_i)),
    
            out_file.write('\nExact Fail\t')
            print '\nExact Fail\t',
            for i in range(num_iter):
                out_file.write('%d\t' % (exact_f[i]+exact_f_i[i]))
                print '%d\t' % (exact_f[i]+exact_f_i[i]),
            out_file.write('%d\t' % (max(exact_f) + max(exact_f_i)))
            print '%d\t' % (max(exact_f) + max(exact_f_i)),
    
            out_file.write('\nSuccess Fixes\t')
            print '\nSuccess Fixes \t',
            for i in range(num_iter):
                out_file.write('%d\t' % (successes[i] + successes_i[i]))
                print '%d\t' % (successes[i] + successes_i[i]),
            out_file.write('%d\t' % (sum(successes) + sum(successes_i)))
            print '%d\t' % (sum(successes) + sum(successes_i)),
    
            out_file.write('\nFailure Fixes\t')
            print '\nFailure Fixes  \t',
            for i in range(num_iter):
                out_file.write('%d\t' % (failures[i] + failures_i[i]))
                print '%d\t' % (failures[i] + failures_i[i]),
            out_file.write('%d\t' % (sum(failures) + sum(failures_i)))
            print '%d\t' % (sum(failures) + sum(failures_i))
            
        out_file.write('\n\n')
        print '\n'

    writeReads(headers,reads,out_read_file)

def fixKmerIter(kmer_file,headers,reads,k,indel_size,num_iter,output_file_str,out_read_file,indel_state):
    '''Performs iterative error correction without the genome'''
    kmers,indices = getKmers(kmer_file)
    exacts = [0 for x in range(num_iter)]
    changes = [0 for x in range(num_iter)]
    exacts_i = []
    changes_i = []
    
    if indel_state == 1:
        indel_str = "ins"
        delOrIns = False
    elif indel_state == -1:
        indel_str = "del"
        delOrIns = True
    elif indel_state == 0:
        indel_str = 'both'
        delOrIns = True
        
        exacts_i = [0 for x in range(num_iter)]
        changes_i = [0 for x in range(num_iter)]
    
    for i in range(num_iter):
        print 'Iter',i
        for km in kmers:
            reads,ex,ch = fixKmer(km,reads,k,indel_size,delOrIns)
            exacts[i] += ex
            changes[i] += ch
            if indel_state == 0:
                reads,ex,ch = fixKmer(km,reads,k,indel_size,False)
                exacts_i[i] += ex
                changes_i[i] += ch
    writeCorrStats(indel_state,exacts,changes,exacts_i,changes_i,output_file_str)
    writeReads(headers,reads,out_read_file)


def fixKmer(kmer,reads,k,indel_size,delOrIns):
    #print kmer
    if delOrIns:
        curr_kmers = genDelKmers(kmer,indel_size)
        adj_k = k - indel_size
        out_str = 'del'
    else:
        curr_kmers = genInsKmers(kmer,indel_size)
        adj_k = k + indel_size
        out_str = 'ins'
        
    new_reads = reads[:]
    changes = [0 for x in range(indel_size+1)]
    for i in range(len(reads)):
        for j in range(len(reads[i])-adj_k+1):
            found_short = False
            for m in range(indel_size+1):
                if delOrIns:
                    adj_m = m
                else:
                    adj_m = -1*m
                r_kmer = reads[i][j:j+k-adj_m]
                if j < len(reads[i])-k+1+adj_m and r_kmer in curr_kmers[m]:
                    changes[m] += 1
                    #if adj_m == 0:
                    #    print 'exact match at %d of %d' % (j,len(reads[i]))
                    if r_kmer[m:] == kmer:
                        changes[m] -= 1
                        changes[0] += 1
                    else:
                        new_reads[i] = reads[i][:j]+kmer+reads[i][k-adj_m+j:]
                    #    print 'kmer substituted for %s' % (r_kmer)
                    #    print '\t%s=%d, %d of %d' % (out_str,m,j,len(reads[i]))
                    found_short = True
                    break
            if found_short:
                break
    
    return (new_reads,changes[0],sum(changes[1:]))


def getKmers(kmer_file):
    '''Given a file of kmers and indices, return a list of all the kmers and indices'''
    rows = []
    k_file = open(kmer_file,'r')
    for i,line in enumerate(k_file):
        if i > 0:
            rows.append(line.split(' '))
    
    kmers,indices = zip(*rows)
    indices = map(int,indices)
    return (kmers,indices)


def writeCorrStats(indel_state,exacts,changes,exacts_i,changes_i,output_file_str):
    out_file = open(output_file_str,'w')
    out_file.write('Summary Statistics:\n')
    out_file.write('Iteration\t')
    for i in range(num_iter):
        out_file.write('%d\t' % (i+1))
    out_file.write('Tot\t')

    if indel_state == -1 or indel_state == 0:
        out_file.write('\nDeletions:')
    elif indel_state == 1:
        out_file.write('\nInsertions:')
    
    out_file.write('\nExacts\t')
    for i in range(num_iter):
        out_file.write('%d\t' % exacts[i])
    out_file.write('%d\t' % max(exacts))

    out_file.write('\nChanges\t')
    for i in range(num_iter):
        out_file.write('%d\t' % changes[i])
    out_file.write('%d\t' % sum(changes))
    
    if indel_state == 0:
        out_file.write('\nInsertions:')
        
        out_file.write('\nExacts\t')
        for i in range(num_iter):
            out_file.write('%d\t' % exacts_i[i])
        out_file.write('%d\t' % max(exacts_i))

        out_file.write('\nChanges\t')
        for i in range(num_iter):
            out_file.write('%d\t' % changes_i[i])
        out_file.write('%d\t' % sum(changes_i))
        
        out_file.write('\nTotal:')
        
        out_file.write('\nExacts\t')
        for i in range(num_iter):
            out_file.write('%d\t' % (exacts[i]+exacts_i[i]))
        out_file.write('%d\t' % (max(exacts) + max(exacts_i)))

        out_file.write('\nChanges\t')
        for i in range(num_iter):
            out_file.write('%d\t' % (changes[i] + changes_i[i]))
        out_file.write('%d\t' % (sum(changes) + sum(changes_i)))
        
    out_file.write('\n\n')
    out_file.close()

def delFix(reads,genome,k,del_size,offset,kmer,index,out_file_str):
    '''Similar to delCheck but rather than printing and writing the resulting
    kmers with deletions, it overwrites them with the original kmer - error
    correcting those deletions, and it returns the new set of reads with 
    corrected kmers. kmer is kmer found at index, not necessarily centered'''
    out_file = open(out_file_str,'a')
    g_mid = int(len(genome)/2)
    kmer = genome[index:index+k]
    out_file.write('kmer=%s, i=%d\n' % (kmer,index))
    print 'kmer=%s, i=%d' % (kmer,index)
    
    curr_kmers = genDelKmers(kmer,del_size)
    
    #read_changes = []
    new_reads = reads[:]
    success = [0 for x in range(del_size+1)]
    failure = [0 for x in range(del_size+1)]
    for i in range(len(reads)):
        for j in range(len(reads[i])-k+1+del_size):
            found_short = False
            for m in range(del_size+1):
                r_kmer = reads[i][j:j+k-m]
                if j < len(reads[i])-k+1+m and r_kmer in curr_kmers[m]:
                    if abs(index - j) <= abs(len(reads[i])-len(genome))+offset:
                        success[m] += 1
                        if m == 0:
                            out_file.write('exact match at %d of %d\n' % (j,len(reads[i])))
                            print 'exact match at %d of %d' % (j,len(reads[i]))
                        else:
                            #read_changes.append([i,j,m])
                            new_reads[i] = reads[i][:j]+kmer+reads[i][k-m+j:]
                            out_file.write('kmer substituted for %s\n' % (r_kmer))
                            out_file.write('\tdel=%d, %d of %d\n' % (m,j,len(reads[i])))
                            print 'kmer substituted for %s' % (r_kmer)
                            print '\tdel=%d, %d of %d' % (m,j,len(reads[i]))
                    else:
                        failure[m] += 1   
                        out_file.write('%s failed to substitute\n' % (r_kmer))
                        out_file.write('\tdel=%d, %d of %d\n' % (m,j,len(reads[i])))
                        print '%s failed to substitute' % (r_kmer)
                        print '\tdel=%d, %d of %d' % (m,j,len(reads[i]))
                    found_short = True
                    break
            if found_short:
                break

    out_file.write('Successes: ')
    print 'Successes: ',
    for i in range(len(success)):
        if success[i] != 0:
            out_file.write('k=%d:%d ' % (i,success[i]))
            print 'k=%d:%d ' % (i,success[i]),
    out_file.write('Tot:%d\n' % sum(success))
    print 'Tot:%d' % sum(success)

    out_file.write('Failures: ')
    print 'Failures: ',
    for i in range(len(failure)):
        if failure[i] != 0:
            out_file.write('k=%d:%d ' % (i,failure[i]))
            print 'k=%d:%d ' % (i,failure[i]),
    out_file.write('Tot:%d\n\n' % sum(failure))
    print 'Tot:%d\n' % sum(failure)
    
    out_file.close()
    #return (read_changes,success[0],failure[0],sum(success[1:]),sum(failure[1:]))
    return (new_reads,success[0],failure[0],sum(success[1:]),sum(failure[1:]))

def indelFix(reads,genome,k,indel_size,offset,kmer,index,out_file_str,delOrIns,fixIncorrect):
    '''Similar to indelCheck but rather than printing and writing the resulting
    kmers with insertions and deletions, it overwrites them with the original kmer - error
    correcting those deletions, and it returns the new set of reads with 
    corrected kmers. kmer is kmer found at index, not necessarily centered. If
    delOrIns is True, then fix deletions. If false, fix insertions. If fixIncorrect
    is true, fix kmers found at the wrong positions.'''
    out_file = open(out_file_str,'a')
    kmer = genome[index:index+k]
    out_file.write('kmer=%s, i=%d\n' % (kmer,index))
    print 'kmer=%s, i=%d' % (kmer,index)
    
    if delOrIns:
        curr_kmers = genDelKmers(kmer,indel_size)
        adj_k = k - indel_size
        out_str = 'del'
    else:
        curr_kmers = genInsKmers(kmer,indel_size)
        adj_k = k + indel_size
        out_str = 'ins'
    
    new_reads = reads[:]
    success = [0 for x in range(indel_size+1)]
    failure = [0 for x in range(indel_size+1)]
    for i in range(len(reads)):
        for j in range(len(reads[i])-adj_k+1):
            found_short = False
            for m in range(indel_size+1):
                if delOrIns:
                    adj_m = m
                else:
                    adj_m = -1*m
                r_kmer = reads[i][j:j+k-adj_m]
                if j < len(reads[i])-k+1+adj_m and r_kmer in curr_kmers[m]:
                    if abs(index - j) <= abs(len(reads[i])-len(genome))+offset:
                        success[m] += 1
                        if adj_m == 0:
                            out_file.write('exact match at %d of %d\n' % (j,len(reads[i])))
                            print 'exact match at %d of %d' % (j,len(reads[i]))
                        elif r_kmer[m:] == kmer:
                            out_file.write('exact match at %d of %d\n' % (j+1,len(reads[i])))
                            print 'exact match at %d of %d' % (j+1,len(reads[i]))
                            success[m] -= 1
                            success[0] += 1
                        else:
                            new_reads[i] = reads[i][:j]+kmer+reads[i][k-adj_m+j:]
                            out_file.write('kmer substituted for %s\n' % (r_kmer))
                            out_file.write('\t%s=%d, %d of %d\n' % (out_str,m,j,len(reads[i])))
                            print 'kmer substituted for %s' % (r_kmer)
                            print '\t%s=%d, %d of %d' % (out_str,m,j,len(reads[i]))
                    else:
                        failure[m] += 1   
                        if fixIncorrect:
                            new_reads[i] = reads[i][:j]+kmer+reads[i][k-adj_m+j:]
                            out_file.write('%s failed but substituted\n' % (r_kmer))
                            out_file.write('\t%s=%d, %d of %d\n' % (out_str,m,j,len(reads[i])))
                            print '%s failed but substituted' % (r_kmer)
                            print '\t%s=%d, %d of %d' % (out_str,m,j,len(reads[i]))
                        else:                        
                            out_file.write('%s failed to substitute\n' % (r_kmer))
                            out_file.write('\t%s=%d, %d of %d\n' % (out_str,m,j,len(reads[i])))
                            print '%s failed to substitute' % (r_kmer)
                            print '\t%s=%d, %d of %d' % (out_str,m,j,len(reads[i]))
                    found_short = True
                    break
            if found_short:
                break

    out_file.write('Successes: ')
    print 'Successes: ',
    for i in range(len(success)):
        if success[i] != 0:
            out_file.write('%s=%d:%d ' % (out_str,i,success[i]))
            print '%s=%d:%d ' % (out_str,i,success[i]),
    out_file.write('Tot:%d\n' % sum(success))
    print 'Tot:%d' % sum(success)

    out_file.write('Failures: ')
    print 'Failures: ',
    for i in range(len(failure)):
        if failure[i] != 0:
            out_file.write('%s=%d:%d ' % (out_str,i,failure[i]))
            print '%s=%d:%d ' % (out_str,i,failure[i]),
    out_file.write('Tot:%d\n\n' % sum(failure))
    print 'Tot:%d\n' % sum(failure)
    
    out_file.close()
    #return (read_changes,success[0],failure[0],sum(success[1:]),sum(failure[1:]))
    return (new_reads,success[0],failure[0],sum(success[1:]),sum(failure[1:]))

def subFix(reads,genome,k,del_size,offset,kmer,index,out_file_str):
    pass

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
            out_file.write(">%s\n%s\n" % (headers[i],reads[i]))

def testKmer(genome,kmer):
    for i in range(len(genome)):
        if genome[i:i+len(kmer)] == kmer:
            print i

if __name__ == "__main__":
    start_time = time.time()
    main()
    print 'Total time',time.time()-start_time
    
def main():
    center = int(sys.argv[1])
    size = int(sys.argv[2])
    min_overlap = size - 100
    read_file = '/home/mshen/research/data/PacBioCLR/PacBio_10kb_CLR_mapped_whole_removed_homopolymer.fasta'
    genome_file = '/home/mshen/research/data/e_coli_genome.fasta'

    extract(center, size, min_overlap)
    return

def extract(center, size, min_overlap, read_file, genome_file):
    m5_file = '/home/mshen/research/data/blasr_all_m5.out'
    read_out = '/home/mshen/research/extracts/extracted_reads_c%d_s%d.fasta' % (center,size)
    genome_out = '/home/mshen/research/extracts/extracted_genome_c%d_s%d.fasta' % (center,size)
    write_out = True
    
    #Extract the genome and reads mapping to the genome at a center and size
    read_dict = getGenomicRegion(m5_file,size,center,min_overlap)
    headers,reads = extractReads(read_dict,read_file,read_out,write_out)
    genome = extractGenome(fold+genome_file,size,center,fold+genome_out,write_out)
    
    return read_out, genome_out