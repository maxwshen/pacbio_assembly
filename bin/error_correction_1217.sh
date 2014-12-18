#! /bin/csh -f
/home/jeyuan/blasr/alignment/bin/blasr $argv[2] $argv[1] -maxMatch 15 -minMatch 10 -bestn 1 -m 5 -out v0_$argv[1].out;
/home/lin/program/process_alignment v0_$argv[1].out v1_$argv[1];
/home/jeyuan/blasr/alignment/bin/blasr $argv[2] v1_$argv[1] -maxMatch 15 -minMatch 10 -bestn 1 -m 5 -out v1_$argv[1].out;
/home/lin/program/process_alignment v1_$argv[1].out v2_$argv[1];
/home/jeyuan/blasr/alignment/bin/blasr $argv[2] v2_$argv[1] -maxMatch 15 -minMatch 10 -bestn 1 -m 5 -out v2_$argv[1].out;
/home/lin/program/process_alignment v2_$argv[1].out v3_$argv[1];
/home/jeyuan/blasr/alignment/bin/blasr $argv[2] v3_$argv[1] -maxMatch 15 -minMatch 10 -bestn 1 -m 5 -out v3_$argv[1].out;
/home/lin/program/process_alignment v3_$argv[1].out corrected_$argv[1];
rm v0_$argv[1].out;
rm v1_$argv[1].out;
rm v2_$argv[1].out;
rm v3_$argv[1].out;
rm  v1_$argv[1];
rm  v2_$argv[1];
rm  v3_$argv[1];
