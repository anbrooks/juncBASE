#!/bin/bash

libdir=$1
inputDir=$2

python ${libdir}preProcess_getASEventReadCounts_by_chr_step2.py -i $inputDir --by_chr

# for scatter job
chrList=()
for f in *.bed
do
  echo $f
  arr=($(echo $f | tr "_" "\n"))
  chrList=("${chrList[@]} "${arr[1]})
done

dirs=($(ls -d ${inputDir}/*/ | xargs -I dir basename dir))
for sampleDir in ${dirs[@]}
do
  for chr in ${chrList}
    do
      expected_out_file=${sampleDir}_${chr}_intron_exon_junction_counts.txt
      #echo $expected_out_file
      #print the params passing to scatter.sh
      echo ${inputDir}/${sampleDir} ${sampleDir}_${chr} tmp_${chr}_preProcess_getASEventReadCounts_step2.bed
    done
done

# for gather job
#echo