#!/bin/sh

INPUT_PATH="../data"
OUTPUT_PATH="../data"


./seq2expr -s $INPUT_PATH/seqs.fa -e $INPUT_PATH/expr.tab -m $INPUT_PATH/factors.wtmx -f $INPUT_PATH/factor_expr.tab -fo $OUTPUT_PATH/obs_pre.txt -o Direct

