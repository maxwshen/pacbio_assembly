#! /bin/csh -f
/home/jeyuan/blasr/alignment/bin/blasr $argv[2] $argv[1] -maxMatch 15 -minMatch 10 -bestn 1 -m 5 -out v0_$argv[1].out;
/home/lin/program/process_alignment_0108 v0_$argv[1].out corrected_$argv[1];
rm v0_$argv[1].out;
