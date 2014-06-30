#!/bin/sh

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/starting_point"


 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_untrained.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_weak.txt \
-pp $OUTPUT_PATH/par_untrained_weak.par \
-i $INPUT_PATH/factor_info.txt \
-p testparam.save \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

# Example Falgs
# -s ../data/seqs_segal.fa  -e ../data/expr_segal.tab -m ../data/factors_propria.wtmx -c ../data/coop.txt -f ../data/factor_expr_segal.tab -fo ../data/obs_pre.txt -i ../data/factor_info.txt -r ../data/rep.txt -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ 
