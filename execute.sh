#!/bin/sh

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/test"

 echo "Run 0"

./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_christophe.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_pseudocount.txt \
-pp $OUTPUT_PATH/par_pseudocount.par \
-p testparam.save \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

for i in `seq 1 3`
do
 echo "Run $i"
./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_christophe.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_pseudocount.txt \
-pp $OUTPUT_PATH/par_pseudocount.par \
-p param.save \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite
done

# Example Falgs
# -s ../data/Input/seqs_segal.fa  -e ../data/Input/expr_segal.tab -m ../data/Input/factors_christophe.wtmx -c ../data/Input/coop.txt -f ../data/Input/factor_expr_segal.tab -fo ../data/obs_pre.txt -i ../data/Input/factor_info.txt -r ../data/Input/rep.txt -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ valgrind --leak-check=yes -p param.save \-p param.save \

