#!/bin/bash

libdir=$1
input=$2
name=$3
tmp_jcn_file=$4
min_overhang=$5

python ${libdir}/preProcess_getASEventReadCounts_step3.py -i ${input} -n ${name} -t ${tmp_jcn_file} --min_overhang ${min_overhang}