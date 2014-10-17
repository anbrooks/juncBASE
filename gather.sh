#!/bin/bash
. /broad/tools/scripts/useuse
reuse -q Python-2.7

libdir=$1
samples=$2
jcn_seq_len=$3
sample_set_name=$4
thresh=$5
delta_thresh=$6
mt_correction=$7

set -e

chrList=()
for f in *.bed
do
  arr=($(echo $f | tr "_" "\n"))
  chrList=("${chrList[@]} "${arr[1]})
done

#createAS_CountTables.py
echo "running python ${libdir}createAS_CountTables.py -d $(pwd) -s ${samples} -o $(pwd)/tmp_createAS_CountTables_${chr} --which_chr ${chr} --jcn_seq_len ${jcn_seq_len}"
mkdir -p countTables
for chr in ${chrList}
do
  python ${libdir}createAS_CountTables.py -d $(pwd) -s ${samples} -o countTables/tmp_createAS_CountTables_${chr} --which_chr ${chr} --jcn_seq_len ${jcn_seq_len}
done

python ${libdir}combine_createAS_CountTables_by_chr.py -d countTables/ -o ${sample_set_name}

python ${libdir}makeCleanIR_JuncBASE.py --in_prefix ${sample_set_name}

python ${libdir}compareSampleSets.py --in_prefix ${sample_set_name}_IR_cleaned --all_psi_output ${sample_set_name}_IR_cleaned_allAS.txt --thresh ${thresh} --delta_thresh ${delta_thresh} --as_only --mt_correction ${mt_correction}