#!/bin/bash

libdir=$1
inputDir=$2
min_overhang=$3
txt_db1=$4
txt_db2=$5
txt_db3=$6
jcn_seq_len=$7
sqlite_db_dir=$8
sample_set_name=$9
thresh=${10}
delta_thresh=${11}
mt_correction=${12}

python ${libdir}preProcess_getASEventReadCounts_by_chr_step2.py -i $inputDir --by_chr

# for scatter job
chrList=()
for f in *.bed
do
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
      echo $libdir ${inputDir}/${sampleDir} ${sampleDir} ${chr} ${inputDir}/tmp_${chr}_preProcess_getASEventReadCounts_step2.bed ${min_overhang} ${inputDir}/${sampleDir}/${sampleDir}_${chr}/${sampleDir}_${chr} ${txt_db1} ${txt_db2} ${txt_db3} ${jcn_seq_len} ${sqlite_db_dir}
    done
done

# for gather job
printf '%s\n' "${dirs[@]}" > $(pwd)/samples.txt
echo $libdir $(pwd)/samples.txt ${jcn_seq_len} ${sample_set_name} ${thresh} ${delta_thresh} ${mt_correction}
