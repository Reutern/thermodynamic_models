#!/bin/sh

INPUT_PATH="../data"
OUTPUT_PATH="../data"


./seq2expr \
-s $INPUT_PATH/seqs_segal.fa  \
-e $INPUT_PATH/expr_segal.tab \
-m $INPUT_PATH/factors_propria.wtmx \
-c $INPUT_PATH/coop.txt \
-f $INPUT_PATH/factor_expr_segal.tab \
-fo $OUTPUT_PATH/obs_pre.txt \
-i $OUTPUT_PATH/factor_info.txt \
-r $OUTPUT_PATH/rep.txt \
-o Direct	# modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limite


