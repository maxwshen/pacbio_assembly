# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:44:52 2014

@author: Jeffrey
"""

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
            for i in range(k-j):
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
            for i in range(k-j):
                for m in range(len(alpha)):
                    next_kmers.append(curr[:i]+alpha[m]+curr[i:])
        ins_kmers[j+1] = set(next_kmers[:])
        next_kmers = []
    return ins_kmers