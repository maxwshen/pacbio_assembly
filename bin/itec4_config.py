# A parameter config file for itec4.py

prior = '/home/yu/max/research/'
contigs_fold = prior + '/contigs_20kb_full_18/'
overlap_accuracy_cutoff = 75    # .
overlap_length_cutoff = 7000     # .
overlap_accuracy_cutoff_consensus = 98
overlap_length_cutoff_consensus = 7000
num_attempts = 1                # Number of times to try nhood extension.
support_cutoff = 70             # CANDIDATE: Required pct accuracy for support to count
support_ratio = 0.6             # CANDIDATE: Required support for a chosen read from other candidates
limit_km_times_total = 5        # How many times to attempt k-mer matching extension per direction
km_k = 15                       # .
km_cutoff = 100                 # .
support_dist_cutoff = 100000    # CONSENSUS: Bp. length, acceptable support distance from end of consensus
support_t = 3                   # CONSENSUS: Req. # reads to support a position to determine farthest support
nhood_header_limit = float('inf')         # .
nhood_it_limit = 3              # .
n21ratio_cutoff = 0.05             # If n2/n1 is greater than this, find another consensus
blasr_exe = 'blasr'
blasr_zero = 4      # 0 on debruijn, 4 on Yu's computer
blasr_zero_len = 8  # 0 on debruijn, 8 on Yu's computer
blasr_options = '-bestn 1 -m 1'   # Concise output
e_coli_genome = '/home/yu/data/ecoli_consensus_mark.fasta'
ec_prefix = 'C0419_'
use_ecs = False