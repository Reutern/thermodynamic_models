#!/bin/sh

red='\033[0;31m'
yellow='\033[0;33m'
NC='\033[0m' # No Color

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/pwm_comparison"
 

echo " "
 echo "${yellow} untrained Run 0${NC}"
 echo " "

./seq2expr \
-s $INPUT_PATH/seqs.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_untrained_pc.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_full_cic.tab \
-fo $OUTPUT_PATH/obs_pre_untrained.txt \
-pp $OUTPUT_PATH/par_untrained.par \
-o Direct 		# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite


echo " "
 echo "${yellow} wolfe Run 0${NC}"
 echo " "

./seq2expr \
-s $INPUT_PATH/seqs.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_wolfe.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_full_cic.tab \
-fo $OUTPUT_PATH/obs_pre_wolfe.txt \
-pp $OUTPUT_PATH/par_wolfe.par \
-o Direct 		# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite


for i in `seq 1 0`
do
 echo "${yellow} fa based Run $i ${NC}"
 echo " "
./seq2expr \
-s $INPUT_PATH/seqs.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_fa.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_full_cic.tab \
-fo $OUTPUT_PATH/obs_pre_fa.txt \
-pp $OUTPUT_PATH/par_fa.par \
-p param.save \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite
done


# Example Falgs
# -s ../data/Input/seqs.fa  -e ../data/Input/expr_smooth.tab -m ../data/Input/factors_trained.wtmx -c ../data/Input/coop.txt -f ../data/Input/factor_expr_full_cic.tab -fo ../data/obs_pre.txt -i ../data/Input/factor_info.txt -r ../data/Input/rep.txt -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ valgrind --leak-check=yes -p param.save \-p param.save \ -p $OUTPUT_PATH/par_trained_06.par \ 
#-oc $OUTPUT_PATH/occ_pre_weak_ut_03_07.txt \ -ff $INPUT_PATH/free_fix_indicator_file.txt \

