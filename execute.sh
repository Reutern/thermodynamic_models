#!/bin/sh

red='\033[0;31m'
yellow='\033[0;33m'
NC='\033[0m' # No Color

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/old_pwm"
 

echo " "
 echo "${yellow} new cad test${NC}"
 echo " "

./seq2expr \
-s $INPUT_PATH/seqs.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_fa_old.wtmx \
-c $INPUT_PATH/coop.txt \
-p param.save \
-f $INPUT_PATH/factor_expr_full_cic.tab \
-fo $OUTPUT_PATH/obs_pre_test_3.txt \
-pp $OUTPUT_PATH/par_test_3.par \
-o Direct 		# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite



# Example Falgs
# -s ../data/Input/seqs.fa  -e ../data/Input/expr_smooth.tab -m ../data/Input/factors_trained.wtmx -c ../data/Input/coop.txt -f ../data/Input/factor_expr_full_cic.tab -fo ../data/obs_pre.txt -i ../data/Input/factor_info.txt -r ../data/Input/rep.txt -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ valgrind --leak-check=yes -p param.save \-p param.save \ -p $OUTPUT_PATH/par_trained_06.par \ 
#-oc $OUTPUT_PATH/occ_pre_weak_ut_03_07.txt \ -ff $INPUT_PATH/free_fix_indicator_file.txt \

