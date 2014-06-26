#!/bin/sh

INPUT_PATH="../data/Input"
OUTPUT_PATH="../data/SSE_06_ch_deletion"

echo "bcd"

 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_bcd.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_bcd.txt \
-pp $OUTPUT_PATH/par_ch_bcd.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "cad"


 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_cad.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_cad.txt \
-pp $OUTPUT_PATH/par_ch_cad.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "gt"


 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_gt.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_gt.txt \
-pp $OUTPUT_PATH/par_ch_gt.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "kr"

 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_kr.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_kr.txt \
-pp $OUTPUT_PATH/par_ch_kr.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "kni"

 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_kni.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_kni.txt \
-pp $OUTPUT_PATH/par_ch_kni.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "tll"

 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_tll.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_tll.txt \
-pp $OUTPUT_PATH/par_ch_tll.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

echo "hb"

 ./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_smooth.tab  \
-m $INPUT_PATH/factors_ch_hb.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre_ch_hb.txt \
-pp $OUTPUT_PATH/par_ch_hb.par \
-i $INPUT_PATH/factor_info.txt \
-r $INPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite

# Example Falgs
# -s ../data/seqs_segal.fa  -e ../data/expr_segal.tab -m ../data/factors_propria.wtmx -c ../data/coop.txt -f ../data/factor_expr_segal.tab -fo ../data/obs_pre.txt -i ../data/factor_info.txt -r ../data/rep.txt -o Direct
 #-ff $INPUT_PATH/free_fix_indicator_file.txt -p testparam.save \-p testparam.save \ 
