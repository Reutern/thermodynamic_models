#!/bin/sh

red='\033[0;31m'
yellow='\033[0;33m'
NC='\033[0m' # No Color

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/extended"
 

echo " "
 echo "${yellow} Extended data set${NC}"
 echo " "

./seq2expr \
-s $INPUT_PATH/seqs_extended.fa  \
-e $INPUT_PATH/expr_extended.tab  \
-m $INPUT_PATH/factors_fa_full.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_full.tab \
-oc $OUTPUT_PATH/occ_pre_scontrol.txt \
-fo $OUTPUT_PATH/obs_pre_control.txt \
-pp $OUTPUT_PATH/par_control.par \
-o Direct 		# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite



# Example Falgs
# -s ../data/Input/seqs_extended.fa  -e ../data/Input/expr_extended.tab -m ../data/Input/factors_fa.wtmx -c ../data/Input/nocoop.txt -f ../data/Input/factor_expr_old.tab -fo ../data/obs_pre.txt  -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ valgrind --leak-check=yes -p param.save \-p param.save \ -p $OUTPUT_PATH/par_trained_06.par \ 
#-oc $OUTPUT_PATH/occ_pre_weak_ut_03_07.txt \ -ff $INPUT_PATH/free_fix_indicator_file.txt \

