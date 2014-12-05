#! /bin/csh -f
/home/lin/max_data_test/preprocess $argv[1];
mkdir dir_$argv[1]; 
/home/jeyuan/blasr/alignment/bin/blasr v0_$argv[1] v0_$argv[1] -maxMatch 15 -minMatch 10 -bestn 50 -insertion 1 -deletion 1 -match -10 -m 5 -out v0_$argv[1].out;
sed -i 's/\// /' v0_$argv[1].out;
/home/lin/max_data_test/process_alignment_strict_v0 v0_$argv[1].out dir_$argv[1];
cat dir_$argv[1]/*_v1.fasta > v1_$argv[1];
/home/jeyuan/blasr/alignment/bin/blasr v1_$argv[1] v1_$argv[1] -maxMatch 15 -minMatch 10 -bestn 50 -insertion 1 -deletion 1 -match -10 -m 5 -out v1_$argv[1].out;
sed -i 's/\// /' v1_$argv[1].out;
/home/lin/max_data_test/process_alignment_strict_v1 v1_$argv[1].out dir_$argv[1];
cat dir_$argv[1]/*_v2.fasta > v2_$argv[1];
/home/jeyuan/blasr/alignment/bin/blasr v2_$argv[1] v2_$argv[1] -maxMatch 15 -minMatch 10 -bestn 50 -insertion 1 -deletion 1 -match -10 -m 5 -out v2_$argv[1].out;
sed -i 's/\// /' v2_$argv[1].out;
/home/lin/max_data_test/process_alignment_strict_v2 v2_$argv[1].out corrected_$argv[1];
rm -r dir_$argv[1];
rm v0_$argv[1]*;
rm v1_$argv[1]*;
rm v2_$argv[1]*;
