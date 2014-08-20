# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 15:03:29 2014

@author: Jeffrey
"""

import random,time,copy
import numpy as np
from subprocess import call
#import matplotlib.pyplot as plt
import sys,os,getopt

def simulateGenome(N,alpha):
    '''Given a list containing the alphabet of possible letters,
    generates and returns a single sequence of length N'''
    header = 'simulated_genome/0_%d' % N
    return header,''.join(random.choice(alpha) for x in xrange(N))

def simulateReads(genome,cov,read_len,err_mat,alpha,start_time,no_shift):
    '''Given a genome, a coverage, uniform read length, and a list
    of error rates summing to 1 in the order
    ['Match','Mismatch','Insertion','Deletion'],
    output a set of reads and headers matching these parameters.'''
    N = len(genome)
    num_reads = int(N*cov/read_len)
    headers = []
    header_base = 'read_%d%s/genome_%d_%d/mat%d_mis%d_ins%d_del%d/0_%d'
    reads = []
    total_errs = [0 for x in range(len(err_mat))]
    total_conserv = [0 for x in range(N)]
    for i in xrange(num_reads):
        print 'Read %d,%.3f' % (i,time.time()-start_time)
        if no_shift:
            start = 0
        else:
            start = int(i*N/num_reads)
        end = start+read_len
        if end <= N:
            read,errs,conserv = simulateReadErrs(genome[start:end],err_mat,alpha)
            headers.append(header_base % (i,'',start,end, \
                    errs[0],errs[1],errs[2],errs[3],len(read)))
            reads.append(read)
            for j in range(len(errs)):
                total_errs[j] += errs[j]
            for j in range(start,end):
                total_conserv[j] += conserv[j-start]
        else:
            #Read 1 to the end of genome
            read,errs,conserv = simulateReadErrs(genome[start:N],err_mat,alpha)
            headers.append(header_base % (i,'a',start,N, \
                    errs[0],errs[1],errs[2],errs[3],len(read)))
            reads.append(read)
            for j in range(len(errs)):
                total_errs[j] += errs[j]
            for j in range(start,N):
                total_conserv[j] += conserv[j-start]
            #Read 2 from the start of genome
            read,errs,conserv = simulateReadErrs(genome[0:end-N],err_mat,alpha)
            headers.append(header_base % (i,'b',0,end-N, \
                    errs[0],errs[1],errs[2],errs[3],len(read)))
            reads.append(read)
            for j in range(len(errs)):
                total_errs[j] += errs[j]
            for j in range(0,end-N):
                total_conserv[j] += conserv[j-0]
    
    return headers,reads,total_errs,total_conserv

def simulateReadErrs(read,err_mat,alpha):
    '''Given a read and a list of error rates in the order
    ['Match','Mismatch','Insertion','Deletion'], output a new read with
    these errors randomly incorporated at these rates.'''
    new_read = ''
    errs = [0 for x in range(4)]
    i = 0
    conserv = [0 for x in range(len(read))]
    while i < len(read):
        err_i = weightedProb(err_mat)
        if err_i == 0:
            new_read += read[i]
            conserv[i] = 1
            i += 1
        elif err_i == 1:
            new_read += random.choice(alpha)
            i += 1
        elif err_i == 2:
            new_read += random.choice(alpha)
        elif err_i == 3:
            i += 1
        else:
            raise Exception('Index not in range(4)')
        errs[err_i] += 1
    return new_read,errs,conserv

def weightedProb(probs):
    '''Given a list of probabilities, return the index of the
    one that was sampled randomly from a uniform distribution'''
    cutoffs = np.cumsum(probs)
    return cutoffs.searchsorted(np.random.uniform(0,cutoffs[-1]))

def calculateConserv(m5_file,c_file):
    outputs = readBlasr_m5(m5_file)
    conserv = [0 for x in range(outputs[0][7])]
    for i in range(len(outputs)):
        #8,9,10 = tStart,tEnd,tStrand
        #17,18,19 = qSeq,matchPatt,tSeq
        if outputs[i][10] == '+':
            for j in range(outputs[i][8],outputs[i][9]):
                if outputs[i][19][j-outputs[i][8]] != '-' and outputs[i][18][j-outputs[i][8]] == '|':
                    conserv[j] += 1
        else:
            print 'Read',i,'Negative Strand'
    
    make_hist = True
    writeConserv(conserv,c_file,make_hist)

def main(N,alpha,cov,read_len,err_mat,g_file,r_file,s_file,c_file,start_time,no_shift):
    g_head,genome = simulateGenome(N,alpha)
    print 'Made genome,%.3f' % (time.time()-start_time)
    headers,reads,total_errs,total_conserv = simulateReads(genome,cov,read_len,err_mat,alpha,start_time,no_shift)
    print 'Made reads,%.3f' % (time.time()-start_time)
    writeGenome(g_head,genome,g_file)
    writeReads(headers,reads,r_file)
    writeStats(genome,reads,total_errs,s_file)
    make_hist = True
    writeConserv(total_conserv,c_file,make_hist)

def writeGenome(header,genome,file_name):
    '''Given a header and a genome, write the header and genome to
    the file in FASTA format'''
    out_file = open(file_name,'w')
    out_file.write('>%s\n%s\n' % (header,genome))
    out_file.close()

def writeReads(headers,reads,file_name):
    '''Given a list of headers, reads, and a file to write to, writes
    the reads to the file in FASTA format'''
    out_file = open(file_name,'w')
    for i in range(len(headers)):
        gt = ''
        if headers[i][0] != '>':
            gt = '>'
        out_file.write("%s%s\n%s\n" % (gt,headers[i],reads[i]))
    out_file.close()

def writeStats(genome,reads,total_errs,stats_file):
    '''Given a genome, a set of reads, the number of errors, and a file,
    write basic stats for this set of reads'''
    out_file = open(stats_file,'w')
    totalReadLen = 0
    for j in range(len(reads)):
        totalReadLen += len(reads[j])

    out_file.write('Statistics:\n\n')
    out_file.write('Length of genome: %d\n' % len(genome))
    out_file.write('Total # of reads: %d\n' % len(reads))
    out_file.write('Average read length: %.3f\n' % (totalReadLen/float(len(reads))))
    out_file.write('Average read coverage: %.3f\n\n' % (totalReadLen/float(len(genome))))
    out_file.write('\tMat\tMis\tIns\tDel\n')
    out_file.write('Tot:')
    for i in range(len(total_errs)):
        out_file.write('\t%d' % total_errs[i])
    out_file.write('\nAvg:')
    for i in range(len(total_errs)):
        out_file.write('\t%.3f' % (total_errs[i]/float(len(reads))))
    out_file.write('\n\n')
    
    out_file.flush()
    out_file.close()

def writeConserv(total_conserv,conserv_file,make_hist):
    '''Write out the conservation score of the reads at each position
    of the genome (conservation refers to the number of reads that
    maintain the right nt at a specific location in the genome)'''
    binRange = range(min(total_conserv),max(total_conserv)+2)
    if make_hist:
        fig = plt.figure(1)
        plt.clf()
        ax = fig.add_subplot(111)        
        ax.hist(total_conserv,binRange)
        plt.show()
    c_file = open(conserv_file,'w')
    freq,bins = np.histogram(total_conserv,binRange)
    c_file.write('Bin\t')
    for i in range(len(bins)-1):
        c_file.write('%d\t' % bins[i])
    c_file.write('\nNum\t')
    for i in range(len(freq)):
        c_file.write('%d\t' % freq[i])
    c_file.write('\n')
    
    
    for i in range(len(total_conserv)):
        if i % 50 == 0:
            c_file.write('\n%d\t' % i)
        fill_zero = ''
        if total_conserv[i] < 10:
            fill_zero = '0'
        c_file.write('%s%d ' % (fill_zero,total_conserv[i]))
    c_file.write('\n\nTotal:%d, Avg:%.3f\n' % (sum(total_conserv),sum(total_conserv)/float(len(total_conserv))))
    c_file.flush()
    c_file.close()

def genAllKmers(k,alpha,out):
    '''Generates all kmers with length k from an alphabet, alpha'''
    if k == 0:
        return ['']
    elif out == []:
        return genAllKmers(k,alpha,alpha)
    else:
        if len(out[0]) == k:
            return out
        else:
            kmers = []
            for s in out:
                for b in alpha:
                    kmers.append(s+b)
            return genAllKmers(k,alpha,kmers)

def localMultiplicity(contig_file,read_file,genome_file,base_k,add_k,out_file,gi,contig_i,num_iter,num_ins):
    '''Consider the supercontigs of high confidence from contig_file,
    consider all possible 15-mers constructed by extending the last
    10 bp of the supercontig by 5, and calculate the multiplicity
    of these 15-mers in the reads. Find the 15-mers of highest multiplicity
    and determine whether they match the correct 15-mers in the genome.'''
    headers,reads = readFASTA(read_file)
    genome = readGenome(genome_file)
    contigs = readLines(contig_file,0,0)
    alpha = ['A','C','T','G']
    allK = genAllKmers(add_k,alpha,[])
    #contigs[0] starts at genomic index 0, contigs[1] ends at end of genome
    prefix = contigs[contig_i[0]][len(contigs[contig_i[0]])-base_k:]
    suffix = contigs[contig_i[1]][:base_k]
    
    pre_gi = gi[0]
    suf_gi = gi[1]
    o_file = open(out_file,'w')
    o_file.close()
    
    tot = 0
    for it in range(num_iter):
        true_pre = genome[pre_gi-base_k:pre_gi+add_k]
        true_suf = genome[suf_gi-add_k:suf_gi+base_k]
        #all_pre = {prefix+x:[0]*(num_ins+1) for x in allK}
        #all_suf = {x+suffix:[0]*(num_ins+1) for x in allK}
        all_pre = {}
        all_suf = {}
        for x in allK:
            all_pre[prefix+x] = [0]*(num_ins+1)
            all_suf[x+suffix] = [0]*(num_ins+1)
        
        
        pre_hammer = {}
        for p in all_pre:
            '''for ins in genInsKmers(p,1,True)[1]:
                if ins not in pre_hammer:
                    pre_hammer[ins] = [p]
                else:
                    pre_hammer[ins] += [p]'''
            
            pre_ins_kmers = genInsKmers(p,num_ins,True)
            for m in range(num_ins+1):
                for ins in pre_ins_kmers[m]:
                    if ins not in pre_hammer:
                        pre_hammer[ins] = [p]
                    else:
                        pre_hammer[ins] += [p]
        
        suf_hammer = {}
        for s in all_suf:
            '''for ins in genInsKmers(s,1,True)[1]:
                if ins not in suf_hammer:
                    suf_hammer[ins] = [s]
                else:
                    suf_hammer[ins] += [s]'''
            
            suf_ins_kmers = genInsKmers(s,num_ins,True)
            for m in range(num_ins+1):
                for ins in suf_ins_kmers[m]:
                    if ins not in suf_hammer:
                        suf_hammer[ins] = [s]
                    else:
                        suf_hammer[ins] += [s]
        
        k = base_k + add_k
        for i in range(len(reads)):
            for j in range(len(reads[i])-k+1):
                '''kmer = reads[i][j:j+k]
                if kmer in all_pre:
                    all_pre[kmer][0] += 1
                if kmer in all_suf:
                    all_suf[kmer][0] += 1
                if j < len(reads[i])-k:
                    ins_kmer = reads[i][j:j+k+1]
                    if ins_kmer in pre_hammer:
                        for km in pre_hammer[ins_kmer]:
                            all_pre[km][1] += 1
                    if ins_kmer in suf_hammer:
                        for km in suf_hammer[ins_kmer]:
                            all_suf[km][1] += 1'''
                
                for m in range(num_ins+1):
                    if j < len(reads[i])-(k+m)+1:
                        kmer = reads[i][j:j+k+m]
                        if kmer in pre_hammer:
                            for km in pre_hammer[kmer]:
                                all_pre[km][m] += 1
                        if kmer in suf_hammer:
                            for km in suf_hammer[kmer]:
                                all_suf[km][m] += 1
                    
        o_file = open(out_file,'a')
        o_file.write('i = %d\n\n' % it)
        o_file.write('Prefix, %s + all %d-mers; True = %s\n' % (prefix,add_k,true_pre))
        max_pre = ''
        mp = 0
        mp_b = 0
        max_suf = ''
        ms = 0
        ms_b = 0
        for key in sorted(all_pre,key=lambda x: sum(all_pre[x]),reverse=True):
            if not max_pre:
                max_pre = key
                mp = sum(all_pre[key])
                mp_b = all_pre[key][0]
            elif sum(all_pre[key]) == mp and all_pre[key][0] > mp_b:
                print 'tie'
                max_pre = key
                mp = sum(all_pre[key])
                mp_b = all_pre[key][0]
            if sum(all_pre[key]) != 0:
                o_file.write('%s %d\n' % (key,sum(all_pre[key])))
        if max_pre == true_pre:
            tot += 1
            o_file.write('Prefix Success!\n')
            print '%d PS!' % it
            
        o_file.write('\nSuffix, all %d-mers + %s; True = %s\n' % (add_k,suffix,true_suf))
        for key in sorted(all_suf,key=lambda x: sum(all_suf[x]),reverse=True):
            if not max_suf:
                max_suf = key
                ms = sum(all_suf[key])
                ms_b = all_suf[key][0]
            elif sum(all_suf[key]) == ms and all_suf[key][0] > ms_b:
                print 'tie'
                max_suf = key
                ms = sum(all_suf[key])
                ms_b = all_suf[key][0]
            if sum(all_suf[key]) != 0:
                o_file.write('%s %d\n' % (key,sum(all_suf[key])))
        if max_suf == true_suf:
            tot += 1
            o_file.write('Suffix Success!\n')
            print '%d SS!' % it
            
        o_file.write('\n')
        o_file.close()
        
        pre_gi += 5
        suf_gi -= 5
        prefix = true_pre[5:]
        suffix = true_suf[:10]
    print '%d / %d' % (tot,num_iter*2)
        

def genInsKmers(kmer,ins_size,include_end):
    '''Generate a list of all kmers made from ins <= ins_size, where
    ins_kmers[0] = no insertions in k
    ins_kmers[1] = 1 insertion in k ... up to ins_kmers[ins_size]
    if include_end, then insertion will be allowed at end position as well'''
    end = 1
    if include_end:
        end = 0
    k = len(kmer)
    alpha = ['A','C','G','T']
    ins_kmers = [set() for x in range(ins_size+1)]
    ins_kmers[0].add(kmer)
    for j in range(ins_size):
        for curr in ins_kmers[j]:
            for i in range(end,k+j+1-end):
                for m in range(len(alpha)):
                    ins_kmers[j+1].add(curr[:i]+alpha[m]+curr[i:])
    return ins_kmers

def findHighDegreeKmers(reads_file, kmers_file, _k, cutoff):
    _kplus = _k + 1
    #_kminus = _k - 1
    isdna = False
    readcount = 0
    kmers = dict()
    kplusmers = dict()
    #kminusmers = dict()

    f = open(reads_file,'r')
    for i, line in enumerate(f):
        if isdna:
            isdna = False
            dna = line.strip()
            for j in range(len(dna) - _k + 1):
                kmer = dna[j:j+_k]
                if kmer in kmers:
                    kmers[kmer] += 1
                else:
                    kmers[kmer] = 1
            '''for j in range(len(dna) - _kminus + 1):
                kmer = dna[j:j+_kminus]
                if kmer in kminusmers:
                    kminusmers[kmer] += 1
                else:
                    kminusmers[kmer] = 1'''
            for j in range(len(dna) - _kplus + 1):
                kmer = dna[j:j+_kplus]
                if kmer in kplusmers:
                    kplusmers[kmer] += 1
                else:
                    kplusmers[kmer] = 1
        if line[0] == '>' or line[0] == '@':
            readcount += 1
            isdna = True
    f.close()

    degrees = dict()
    for kmer in kmers:
        degree = 0
        '''del_kmers = genDelKmers(kmer, 1, False)[1]
        for del_kmer in del_kmers:
            if del_kmer in kminusmers:
                # degree += 1
                degree += kminusmers[del_kmer]'''

        ins_kmers = genInsKmers(kmer, 1, False)[1]
        for ins_kmer in ins_kmers:
            if ins_kmer in kplusmers:
                # degree += 1 
                degree += kplusmers[ins_kmer]

        degrees[kmer] = degree

    numToOutput = cutoff
    num = copy.copy(numToOutput)
    best = set()
    k_file = open(kmers_file,'w')
    for key in sorted(degrees, key=degrees.get, reverse=True):
        if num == 0:
            break
        k_file.write('%s\n' % key)
        #print key, 'Deg =', degrees[key], 't =', kmers[key]
        best.add(key)
        num -= 1
    k_file.close()
    return best


def contigReadSubstitution(contigs,contig_kmers,reads,k,out_file):
    '''Given a list of contigs with high accuracy based on the
    deBruijn graph, a sorted list of high confidence kmers for each contig,
    a list of reads that the contigs were originally
    made from, and the k used for contig_kmers, this function substitutes
    portions of the contigs into the reads using contig_kmers as anchors
    at the ends and returns the new reads.'''
    new_reads = reads[:]
    o_file = open(out_file,'w')
    read_sub_count = [0 for x in range(len(reads))]
    contig_sub_count = [0 for x in range(len(contigs))]
    for r in range(len(reads)):
        kmer_dict = {}
        for i in range(len(reads[r])-k+1):
            kmer_dict[reads[r][i:i+k]] = i
        toSubInd = []
        toSubSeq = []
        for c in range(len(contigs)):
            matches = [] # each entry will be (contig_index, read_index)
            for j in range(len(contig_kmers[c])):
                if contig_kmers[c][j][0] in kmer_dict:
                    matches.append((contig_kmers[c][j][1],kmer_dict[contig_kmers[c][j][0]]))
            #print 'R',r,'C',c,len(matches)
            if len(matches) > 1:
                f = 0
                l = -1
                tog = True
                failed = False
                while f < len(matches) and (matches[f][0] >= matches[l][0] or matches[f][1] >= matches[l][1]) and not failed:
                    print matches,c,r,contigs[c],reads[r]
                    print 'Failed',matches[f],matches[l]
                    if tog:
                        f += 1
                    else:
                        l -= 1
                    tog = not tog
                    failed = True
                if not failed:
                    readSeq = reads[r][matches[f][1]:matches[l][1]+k]
                    conSeq = contigs[c][matches[f][0]:matches[l][0]+k]
                
                    if conSeq != readSeq:
                        toSubInd.append((matches[f][1],matches[l][1],c))
                        toSubSeq.append(conSeq)
                        read_sub_count[r] += 1
                        contig_sub_count[c] += 1
                    
                        o_file.write('Read %d Contig %d\n' % (r,c))
                        o_file.write('read i %d:%d\n' % (toSubInd[-1][0],toSubInd[-1][1]+k))
                        o_file.write('contig i %d:%d\n' % (matches[f][0],matches[l][0]+k))
                        o_file.write('%s\nsubstituted by\n%s\n\n' % (readSeq,conSeq))
                    else:
                        pass
                        #print 'Iden'
                
        if toSubInd:
            o_file.write('\n')
            o_file.flush()
            toSubInd,toSubSeq = zip(*sorted(zip(toSubInd,toSubSeq)))
            curr_i_old = 0
            new_reads[r] = ''
            for j in range(len(toSubInd)):
                start_i,end_i,c = toSubInd[j]
                new_reads[r] += reads[r][curr_i_old:start_i]+toSubSeq[j]
                curr_i_old = end_i
            new_reads[r] += reads[r][curr_i_old:]
    o_file.close()
    return new_reads,read_sub_count,contig_sub_count

def filterContigs(contigs,kmer_dict,k,min_t,min_contig_len):
    '''Given a set of contigs, a dictionary of the degree/t of all kmers
    in the contigs, and a threshold for t and contig length, filter all
    contigs that are below the minimum length and return a list equal in
    length to the remaining contigs, with each entry being a sorted list of
    kmers present above the threshold for that contig with their indices.'''
    
    new_contigs = []
    contig_kmers = []
    for c in contigs:
        if len(c) > min_contig_len:
            new_contigs.append(c)
            contig_kmers.append([])
            for i in range(len(c)-k+1):
                kmer = c[i:i+k]
                if kmer in kmer_dict:
                    if kmer_dict[kmer][1] >= min_t:
                        contig_kmers[-1].append((kmer,i))
                '''else:
                    print "Missing kmer in contig"
                    print kmer,c'''
        else:
            print 'Filtered out contig of length %d' % len(c)
    return new_contigs,contig_kmers


def readBlasr_m5(m5_file_str):
    '''Reads and returns all outputs of a blasr file'''

    readOut = []    
    m5_file = open(m5_file_str,'r')
    #           0   1       2       3   4       5       6   7       8       9   10      11      12      13          14      15  16      17          18          19
    #Output: qName qLength qStart qEnd qStrand space tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
    num_i = [1,2,3,7,8,9,11,12,13,14,15,16]
    for i,line in enumerate(m5_file):
        parts = line.strip().split(" ")
        for j in num_i:
            parts[j] = int(parts[j])
        readOut.append(parts[:])
    m5_file.close()
    return readOut

def getIndelsFromAlign(qSeq,alignSeq,tSeq):
    '''Given an alignment set of 3 sequences, count the total number of
    insertions and deletions (where "ACA---AT" counts as 1 deletion) as well
    as the number of each consecutive indel and return them'''
    
    #format for cons_counts = [Mat,Mis,Ins,Del], each a dict from len -> num
    cons_counts = [{} for x in range(4)]
    curr_count = 1
    prev_i = -1
    curr_i = 0
    for i in range(len(qSeq)):
        if qSeq[i] == '-':
            curr_i = 3
        elif tSeq[i] == '-':
            curr_i = 2
        elif alignSeq[i] == '*':
            curr_i = 1
        elif alignSeq[i] == '|':
            curr_i = 0
        
        if curr_i == prev_i:
            curr_count += 1
        elif prev_i != -1:
            if curr_count not in cons_counts[prev_i]:
                cons_counts[prev_i][curr_count] = 1
            else:
                cons_counts[prev_i][curr_count] += 1
            curr_count = 1
        if i != len(qSeq)-1:
            prev_i = curr_i

    if curr_i == prev_i:
        if curr_count not in cons_counts[prev_i]:
            cons_counts[prev_i][curr_count] = 1
        else:
            cons_counts[prev_i][curr_count] += 1
    
    return cons_counts

def blasrAlignment(query,tfile,in_file,out_file):
    '''Uses blasr to produce the alignments of seq1 and seq2 to the (temp) out_file
    in m 5 format, and then reads the output to find scores. Each query should
    be one read, and the target should be the genome. The output is a list of
    all ['Match','Mis','Ins','Del','NumIns','NumDel','AllErrs','All','Score'] 
    and consecutive indels for the read as28 27 29 30 28 27 30 28 30 29 30 28 30 30 28 28 30 29 27 29 30 29 28 30 30 28 28 30 27 28 30 27 28 27 30 28 29 29 30 30 29 30 28 29 28 30 27 28 27 27 28 29 30 26 28 28 30 29 29 27 30 29 30 28 30 30 28 29 29 28 30 28 28 30 28 30 30 30 28 29 29 28 27 28 29 30 26 28 29 28 27 30 29 29 30 29 30 29 30 26 28 30 28 29 30 30 25 29 30 28 26 30 28 30 29 29 28 30 28 28 28 28 30 28 28 26 29 27 27 26 28 29 29 30 28 29 30 29 29 30 29 28 28 30 28 26 29 30 28 28 30 28 28 30 27 29 30 30 27 30 29 28 30 30 30 30 29 28 29 30 29 28 28 30 29 29 29 29 28 28 28 28 30 28 29 29 28 28 29 27 28 28 28 29 28 27 29 29 29 29 30 29 30 27 28 29 29 28 29 29 28 30 28 29 29 30 27 29 29 29 30 29 29 30 30 29 27 27 29 30 29 28 30 29 30 30 29 28 29 29 29 30 29 30 29 27 27 30 26 30 30 30 28 29 30 27 28 30 28 28 30 30 29 29 26 30 29 29 29 26 30 30 27 30 28 30 29 30 29 29 28 27 29 30 27 30 30 29 28 30 28 30 27 29 28 27 28 30 28 28 28 30 30 26 29 30 29 30 28 29 29 30 30 30 28 29 29 30 29 29 28 26 27 29 28 29 30 26 29 30 26 30 27 28 30 30 28 28 29 28 29 29 29 30 30 29 29 27 26 28 30 28 27 30 27 29 29 30 27 28 28 28 29 29 26 29 27 27 30 27 28 28 29 30 29 29 30 28 30 30 29 29 30 29 30 27 29 30 30 30 30 27 28 28 28 29 29 30 29 28 30 30 30 29 29 28 29 29 28 30 30 29 30 28 29 30 29 28 27 30 27 28 30 30 30 28 29 28 29 29 27 29 29 29 28 28 28 29 29 29 30 29 29 28 27 30 29 29 29 26 28 30 29 29 29 30 28 28 30 27 30 30 28 30 29 28 29 30 27 28 29 29 29 29 28 28 29 30 30 28 29 25 30 29 30 30 30 28 30 30 29 30 29 29 27 29 28 28 28 28 30 29 28 29 30 28 29 30 29 29 29 30 30 27 30 29 30 28 30 29 29 29 30 29 29 27 28 27 28 28 29 28 29 30 29 29 28 30 28 29 29 29 29 30 30 27 30 28 29 30 29 26 28 30 30 30 30 28 30 29 28 28 26 30 27 28 27 30 28 29 29 30 27 29 30 28 30 29 29 28 29 28 28 29 30 29 29 28 28 29 29 30 26 29 29 29 27 28 28 30 28 29 27 30 28 28 29 28 29 28 30 28 29 30 30 29 30 30 30 28 30 30 30 30 29 28 29 30 29 27 29 30 29 28 30 28 28 30 30 29 30 29 29 29 28 28 28 29 29 28 29 27 30 30 30 29 29 28 30 29 25 27 28 29 30 29 30 30 29 27 27 27 26 29 29 27 29 29 29 30 28 28 28 29 29 29 29 28 29 28 29 30 29 29 28 29 29 30 29 29 29 29 30 28 26 28 30 29 29 28 29 28 30 29 28 30 29 29 28 29 30 30 28 30 28 29 29 27 29 30 29 28 29 26 29 28 30 30 29 30 29 29 27 29 27 26 29 27 28 27 29 28 30 28 28 30 28 29 28 30 29 30 30 30 26 30 28 29 29 29 28 30 29 30 29 29 29 29 29 30 30 29 29 28 30 29 29 29 28 25 26 29 29 30 29 27 30 28 29 30 30 29 29 26 29 27 28 30 28 30 29 28 30 26 30 30 29 27 28 30 28 29 29 27 28 29 30 29 29 29 30 29 28 28 30 30 26 30 30 27 28 30 29 30 28 28 26 28 30 30 29 27 28 28 29 29 29 30 26 27 30 30 29 30 29 30 30 30 29 30 29 30 30 30 28 30 29 28 30 28 28 29 27 23 30 30 30 29 27 29 30 29 27 30 29 29 26 29 29 28 30 30 28 29 30 28 29 30 30 28 30 29 29 30 29 27 29 29 27 28 26 29 26 30 30 29 30 27 29 30 29 24 26 28 30 29 30 29 28 28 26 30 28 29 28 28 28 28 28 28 27 28 29 30 30 29 30 27 30 30 26 28 30 29 30 29 30 26 28 30 29 29 28 27 27 29 30 28 30 30 29 29 30 27 30 29 28 28 27 26 28 29 25 29 29 26 30 29 27 30 well as the alignment.
    Uses default blasr parameters.'''
    q_head = 'temp_blasr_read'
    writeReads([q_head],[query],in_file)
    writeBlasr(tfile,in_file,out_file,5)
    #blasr_comm = ['blasr',in_file,tfile,'-bestn','1','-m','5','-out',out_file]
    #call(blasr_comm)
    outputs = readBlasr_m5(out_file)
    if not outputs:
        return ([0 for x in range(9)],[{} for x in range(4)],'','','',0,0,'')
    
    allCounts = [0 for x in range(9)]
    i = 0   #Assumes only 1 read
    for j in range(4):
        allCounts[j] = outputs[i][j+12]
    allCounts[6] = allCounts[1]+allCounts[2]+allCounts[3]
    allCounts[7] = allCounts[6]+allCounts[0]
    allCounts[8] = outputs[i][11]
    qSeq = outputs[i][17]
    alignSeq = outputs[i][18]
    tSeq = outputs[i][19]
    tStart = outputs[i][8]
    tEnd = outputs[i][9]
    tStrand = outputs[i][10]
    #DEAL WITH NEGATIVE STRAND HERE
    
    consCounts = getIndelsFromAlign(qSeq,alignSeq,tSeq)

    for key in consCounts[2]:
        allCounts[4] += consCounts[2][key]
    for key in consCounts[3]:
        allCounts[5] += consCounts[3][key]
    
    return (allCounts,consCounts,qSeq,alignSeq,tSeq,tStart,tEnd,tStrand)

def multBlasrAlign(tInds,aligns,genome,mult_file):
    '''Given the indices for the start, end, and strands for the alignment of all reads
    to the genome for one specific iteration, and the sequences of the alignments
    of these reads to the genome, and the genome itself, write all of the sequences
    aligned to the genome in one multiple alignment.'''
    #Sort both lists by the first index of tInds (tStart)
    tInds,aligns = zip(*sorted(zip(tInds,aligns),key=lambda x: x[0][0]))
    out_lines = ['' for x in range(len(tInds)+1)]
    curr_i = [0 for x in range(len(tInds))]
    i = 0
    while i < len(genome):
        no_ins = True
        for j in range(len(tInds)):
            tStart,tEnd,tStrand = tInds[j]
            qSeq,alignSeq,tSeq = aligns[j]
            if tStart <= curr_i[j] and tEnd > curr_i[j] and tSeq[curr_i[j]-tStart] == '-':
                no_ins = False
                break
        
        if no_ins:
            out_lines[0] += genome[i]
            for j in range(len(tInds)):
                tStart,tEnd,tStrand = tInds[j]
                qSeq,alignSeq,tSeq = aligns[j]
                if tStart > curr_i[j] or tEnd <= curr_i[j]:
                    out_lines[j+1] += ' '
                else:
                    align_i = curr_i[j] - tStart
                    out_lines[j+1] += qSeq[align_i]
                    if tSeq[align_i] != genome[i]:
                        print 'Failed to match','Ind',i,'Read',j,'Align',align_i
                curr_i[j] += 1
            i += 1
            num_ins = 0
        else:
            out_lines[0] += '-'
            for j in range(len(tInds)):
                tStart,tEnd,tStrand = tInds[j]
                qSeq,alignSeq,tSeq = aligns[j]
                if tStart > curr_i[j] or tEnd <= curr_i[j]:
                    out_lines[j+1] += ' '
                else:
                    align_i = curr_i[j] - tStart
                    if tSeq[align_i] == '-':
                        curr_i[j] += 1
                        out_lines[j+1] += qSeq[align_i]
                    else:
                        out_lines[j+1] += '-'
                        
    m_file = open(mult_file,'w')
    for i in range(len(out_lines)):
        m_file.write('%d\t%s\n' % (i,out_lines[i]))
    m_file.flush()
    m_file.close()
    return out_lines


def writeBlasr(genomeFile,inFile,outFile,m,mm=0):
    '''Use blasr to align the reads in inFile with format m'''
    if mm != 0:
        blasr_comm = ['/home/jeyuan/blasr/alignment/bin/blasr',inFile,genomeFile,'-bestn','1','-m',str(m),'-out',outFile,'-minMatch',str(mm)]
    else:
        blasr_comm = ['/home/jeyuan/blasr/alignment/bin/blasr',inFile,genomeFile,'-bestn','1','-m',str(m),'-out',outFile]
    call(blasr_comm)
    
def readLines(line_file,head_size,tail_size):
    '''Simply output all of the lines in a list (with \n stripped), with
    the head and tail excluded.'''
    lines = []
    l_file = open(line_file,'r')
    for i,l in enumerate(l_file):
        if i >= head_size:
            lines.append(l.strip())
    l_file.close()
    while tail_size > 0:
        lines.pop()
        tail_size -= 1
    return lines

def readKmers(kmer_file,head_size,tail_size):
    '''Given a file of kmers, one per line, return a dict of all the kmers.
    header and tail are expected to extend head_size and tail_size lines.'''
    #output format: kmer -> (deg,t,True if correct)
    lines = readLines(kmer_file,head_size,tail_size)
    k_dict = {}
    for l in lines:
        parts = l.split(" ")
        corr = False
        if parts[7] == "correct":
            corr = True
        if parts[0] not in k_dict:
            k_dict[parts[0]] = (int(parts[3]),int(parts[6]),corr)
        else:
            print "Failed: kmer %d showed up multiple times" % parts[0]
    
    return k_dict

def readFASTA(file_name):
    '''Given a FASTA file of reads where each sequence is ON ONE LINE,
    return a list of the headers and a list of the reads of the sequences'''
    seq_file = open(file_name,'r')
    headers = []
    seqs = []
    line = seq_file.readline().strip()
    while line:
        if line[0] == ">":
            headers.append(line)
        else:
            seqs.append(line)
        line = seq_file.readline().strip()
    seq_file.close()
    return (headers,seqs)

def readGenome(file_name):
    '''Given a fasta file with a genome sequence, returns the genome sequence'''
    g_file = open(file_name,'r')
    header = g_file.readline().strip()
    #print header
    genome = ""
    line = g_file.readline().strip()
    while line and line != "":
        genome += line
        line = g_file.readline().strip()
    g_file.close()
    return genome


def writeAlignStats(ac,align_file):
    a_file = open(align_file,'a')
    a_file.write('Match\tMis\tIns\tDel\tNumIns\tNumDel\tAll\tErrRate\tScore\t\n')
    for j in range(len(ac)):
        if j != 7:
            a_file.write('%d\t' % ac[j])
        else:
            a_file.write('%.3f\t' % ac[j])
    a_file.write('\n\n')

def writeConsStats(allCons,cons_file):
    c_file = open(cons_file,'a')
    h = ['Mat','Mis','Ins','Del']
    cons_counts = [{} for x in range(len(allCons[0]))]
    for i in range(len(allCons)):
        for j in range(len(allCons[0])):
            for c in allCons[i][j]:
                if c not in cons_counts[j]:
                    cons_counts[j][c] = allCons[i][j][c]
                else:
                    cons_counts[j][c] += allCons[i][j][c]
    for j in range(len(cons_counts))[::-1]:
        c_file.write('%s\t' % h[j])
        for c in cons_counts[j]:
            c_file.write('%d: %d, ' % (c,cons_counts[j][c]))
        c_file.write('\n')
    c_file.write('\n')
    
def getCounts(align_counts):
    counts = [0 for x in range(len(align_counts[0]))]
    for i in range(len(align_counts)):
        if align_counts[i]:
            for j in range(len(counts)):
                counts[j] += align_counts[i][j]
    errRate = 0.000
    if counts[7] != 0:
        errRate = counts[6] / float(counts[7])
    counts.pop(6)
    counts.insert(7,errRate)
    return counts

def subAlignments(reads,read_file,genome_file,blasr_files,al_stats_file):
    genome = readGenome(genome_file)
    blasr_in,blasr_out,blasr_align_file,mult_file = blasr_files
    tInds = [(0,0,'') for x in range(len(reads))]
    alSeqs = [('','','') for x in range(len(reads))]
    allAligns = [[] for x in range(len(reads))]
    allCons = [[] for x in range(len(reads))]
    
    for i in range(len(reads)):
        ac,cc,readSeq,align,genSeq,tStart,tEnd,tStrand = blasrAlignment(reads[i],genome_file,blasr_in,blasr_out)
        allAligns[i] = ac
        allCons[i] = cc
        tInds[i] = (tStart,tEnd,tStrand)
        alSeqs[i] = (readSeq,align,genSeq)
    
    writeAlignStats(getCounts(allAligns),al_stats_file)
    writeConsStats(allCons,al_stats_file)
    writeBlasr(genome_file,read_file,blasr_align_file,0)
    multBlasrAlign(tInds[:],alSeqs[:],genome,mult_file)

def alignContigs(contigs,genome_file,con_temp_file,con_al_file,al_stats_file):
    headers = ['Contig_%d' % x for x in range(len(contigs))]
    writeReads(headers,contigs,con_temp_file)
    writeBlasr(genome_file,con_temp_file,con_al_file,5,1)
    output = readBlasr_m5(con_al_file)
    al_file = open(al_stats_file,'a')
    for i in range(len(output)):
        al_file.write('Contig %d: %s\n' % (i,output[i][:16]))
    al_file.write('\n')
    al_file.close()

def subContigReads_main(contig_file,kmer_file,read_file,new_read_file,out_file,genome_file,old_blasr_files,new_blasr_files,al_stats_file,con_temp,con_al,k,min_t,min_contig_len):
    k_dict = readKmers(kmer_file,0,0)
    contigs = readLines(contig_file,0,0)
    al_file = open(al_stats_file,'w')
    al_file.close()
    alignContigs(contigs,genome_file,con_temp,con_al,al_stats_file)
    contigs,contig_kmers = filterContigs(contigs,k_dict,k,min_t,min_contig_len)

    headers,reads = readFASTA(read_file)
    al_file = open(al_stats_file,'a')
    al_file.write('Old reads\n')
    al_file.close()
    subAlignments(reads,read_file,genome_file,old_blasr_files,al_stats_file)    

    new_reads,read_counts,contig_counts = contigReadSubstitution(contigs,contig_kmers,reads,k,out_file)
    for i in range(len(headers)):
        headers[i] += '/contig_substituted'
    writeReads(headers,new_reads,new_read_file)
    
    al_file = open(al_stats_file,'a')
    al_file.write('New reads\n')    
    al_file.close()
    subAlignments(new_reads,new_read_file,genome_file,new_blasr_files,al_stats_file)
    
    contig_len = [len(x) for x in contigs]
    
    al_file = open(al_stats_file,'a')
    al_file.write('Read\t')
    for i in range(len(read_counts)):
        al_file.write('%d\t' % i)
    al_file.write('Tot\tAvg\nNum Sub\t')
    for i in range(len(read_counts)):
        al_file.write('%d\t' % read_counts[i])
    al_file.write('%d\t%.3f\n\nContig\t' % (sum(read_counts),sum(read_counts)/float(len(read_counts))))
    for i in range(len(contig_counts)):
        al_file.write('%d\t' % i)
    al_file.write('Tot\tAvg\nLen\t')
    for i in range(len(contig_len)):
        al_file.write('%d\t' % contig_len[i])
    al_file.write('%d\t%.3f\nNum Sub\t' % (sum(contig_len),sum(contig_len)/float(len(contig_len))))
    for i in range(len(contig_counts)):
        al_file.write('%d\t' % contig_counts[i])
    al_file.write('%d\t%.3f\n\n' % (sum(contig_counts),sum(contig_counts)/float(len(contig_counts))))
    al_file.close()

if __name__ == "__main__":
    start_time = time.time()
    N = 1000
    alpha = ['A','C','G','T']
    cov = 30
    read_len = 1000
    # [Mat, Mis, Ins, Del]
    err_mat = [.853,.014,.108,.025]
    g_file = 'sim_nsgenome_1kb.fasta'
    r_file = 'sim_nsreads_1kb.fasta'
    s_file = 'sim_nsstats_1kb.txt'
    c_file = 'sim_nsconserv_1kb.txt'
    no_shift = True
    #main(N,alpha,cov,read_len,err_mat,g_file,r_file,s_file,c_file,start_time,no_shift)
    '''
    gi = 68,887
    contig_i = 9,10
    num_iter = 150
    num_ins = 1
    #contig_file = 'sim_nsregion.s4.t3.11.L1.contigs.txt'
    contig_file = 'sim_nsregion.s0.t3.10.L0.contigs.txt'
    read_file = 'sim_reads_1kb_noshift.fasta'
    genome_file = 'sim_genome_1kb_noshift.fasta'
    local_mult_file = 'sim_noshift_local_mult.txt'
    localMultiplicity(contig_file,read_file,genome_file,10,5,local_mult_file,gi,contig_i,num_iter,num_ins)
    '''
    #m5_file = 'blasr_noshift_1kb_test_m5.out'
    #conserv_file = 'sim_noshift_1kb_conserv.txt'
    #calculateConserv(m5_file,conserv_file)
    
    #contig_file = 'sim_region.15.top1000.contigs.txt'
    #kmer_file = 'sim_region.15.top1000.kmers.txt'
    
    
    fold,contig_file,kmer_file,genome_file,read_file,new_read_file,k,min_contig_len,min_t = ('','','','','','',0,0,0)
    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"hf:c:m:g:r:o:k:l:t:",["fold=","contigs=","kmers=","genome=","inreads=","outreads=","k=","minconlen=","min_t="])
    except getopt.GetoptError:
        print 'GenomeSimulator.py -f <output_fold> -c <contig_file> -m <kmer_file> -g <genome_file> -r <input_reads> -o <output_reads -k <k> [-m <min_contig_len> (10+k)] [-t <min_t_of_kmers> (4)]'
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h' or opt == '--help':
            print 'GenomeSimulator.py -f <output_fold> -c <contig_file> -k <kmer_file> -g <genome_file> -r <input_reads> -o <output_reads -k <k> [-m <min_contig_len> (10+k)] [-t <min_t_of_kmers> (4)]'
            sys.exit()
        elif opt in ('-f','--fold'):
            fold = arg
            if fold and fold[-1] != '/':
                fold += '/'
        elif opt in ('-c','--contigs'):
            contig_file = arg
        elif opt in ('-m','--kmers'):
            kmer_file = arg
        elif opt in ('-g','--genome'):
            genome_file = arg
        elif opt in ('-r','--inreads'):
            read_file = arg
        elif opt in ('-o','--outreads'):
            new_read_file = arg
        elif opt in ('-k','--k'):
            k = int(arg)
        elif opt in ('-m','--mincontig'):
            min_contig_len = int(arg)
        elif opt in ('-t','--min_t'):
            min_t = int(arg)
        
    if min_contig_len == 0:
        min_contig_len = k + 10
    if min_t == 0:
        min_t = 4
    if contig_file == '' or kmer_file == '' or genome_file == '' or read_file == '' or new_read_file == '' or k == 0:
        print 'GenomeSimulator.py -f <output_fold> -c <contig_file> -m <kmer_file> -g <genome_file> -r <input_reads> -o <output_reads -k <k> [-m <min_contig_len> (10+k)] [-t <min_t_of_kmers> (4)]'
        sys.exit(2)
    else:
        if not os.path.exists(fold):
            os.makedirs(fold)

        '''s = 0
        k = 10-s
        # value of t for 14,13,12,...
        t = [3,2,2,3,3,3]
        #c_iters = [(1,2,14),(2,2,13),(3,9,12)]
        #k_iters = [(14,1,2),(13,2,2),(12,3,9)]
        c_iter = (s,t[s],k)
        k_iter = (k,s,t[s])
        contig_file = 'sim_nsregion.s%d.t%d.%d.L0.contigs.txt' % c_iter
        kmer_file = 'sim_nsregion_highdegnodes.%d.s%d.t%d.out' % k_iter
        #k = k_iters[i][0]
        #k = 15
        read_files = ['sim_nsregion_subst_reads_k%d.fasta' % x for x in range(15-s,15+1)[::-1]]
        rd_file_orig = 'sim_reads_1kb_noshift.fasta'
        read_files.insert(0,rd_file_orig)
        read_file = read_files[s]
        genome_file = 'sim_genome_1kb_noshift.fasta'    
        new_read_file = read_files[s+1]'''
        
        out_file = 'sim_nsregion_subst_output_k%d.txt' % k
        blasr_in = 'sim_nsregion_blasr_read.fasta'
        blasr_out = 'sim_nsregion_blasr_m5.out'    
        blasr_align_old = 'sim_nsregion_blasr_k%d_old_m0.out' % k
        blasr_align_new = 'sim_nsregion_blasr_k%d_new_m0.out' % k
        mult_old = 'sim_nsregion_mult_align_k%d_old.out' % k
        mult_new = 'sim_nsregion_mult_align_k%d_new.out' % k
        al_stats_file = 'sim_nsregion_al_k%d.out' % k
        old_blasr_files = (fold+blasr_in,fold+blasr_out,fold+blasr_align_old,fold+mult_old)
        new_blasr_files = (fold+blasr_in,fold+blasr_out,fold+blasr_align_new,fold+mult_new)
        con_temp_file = 'sim_nsregion_contigs_k%d.fasta' % k
        con_al_file = 'sim_nsregion_contig_align_k%d.out' % k
        #min_t = 4
        #min_contig_len = 25
        subContigReads_main(contig_file,kmer_file,read_file, \
                            fold+new_read_file,fold+out_file,genome_file, \
                            old_blasr_files,new_blasr_files,fold+al_stats_file, \
                            fold+con_temp_file,fold+con_al_file,k,min_t,min_contig_len)
        
        print 'Total time',time.time()-start_time