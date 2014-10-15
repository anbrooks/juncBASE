#!/bin/bash
. /broad/tools/scripts/useuse
reuse -q Python-2.7

libdir=$1
samples=$2
chrList=$3
jcn_seq_len=$4
sample_set_name=$5
thresh=$6
delta_thresh=$7
mt_correction=$8

set -e

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