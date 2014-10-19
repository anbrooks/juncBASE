#!/bin/bash

. /broad/tools/scripts/useuse
reuse .biopython-1.61-python-2.7.1-sqlite3-rtrees

set -e

libdir=$1
outputDir=$2
input=$3
sample=$4
chr=$5
tmp_jcn_file=$6
min_overhang=$7
jnc_path_prefix=$8
txt_db1=$9
txt_db2=${10}
txt_db3=${11}
jcn_seq_len=${12}
sqlite_db_dir=${13}

echo "running python ${libdir}preProcess_getASEventReadCounts_step3.py -i ${input} -n ${sample}_${chr} -t ${tmp_jcn_file} --min_overhang ${min_overhang}"
python ${libdir}preProcess_getASEventReadCounts_step3.py -i ${input} -n ${sample}_${chr} -t ${tmp_jcn_file} --min_overhang ${min_overhang}

#the output prefix structure is for the gather step
mkdir -p ${outputDir}/${sample}/${sample}_${chr}
echo "running python ${libdir}/getASEventReadCounts.py --jcn1 ${tmp_jcn_file} --jnc2 ${jnc_path_prefix}_junctions.bed --genome_reads1 ${jnc_path_prefix}_genome_reads.txt.gz
 --genome_reads2 ${jnc_path_prefix}_genome_reads.txt.gz --ie1 ${jnc_path_prefix}_intron_exon_junction_counts.txt --ie2 ${jnc_path_prefix}_intron_exon_junction_counts.txt
 -p ${outputDir}/${sample}/${sample}_${chr}/${sample}_${chr} --txt_db1 ${txt_db1} --txt_db2 ${txt_db2} --txt_db3 ${txt_db3} --jcn_seq_len ${jcn_seq_len} --sqlite_db_dir ${sqlite_db_dir}"
python ${libdir}getASEventReadCounts.py --jcn1 ${tmp_jcn_file} --jcn2 ${jnc_path_prefix}_junctions.bed --genome_reads1 ${jnc_path_prefix}_genome_reads.txt.gz \
 --genome_reads2 ${jnc_path_prefix}_genome_reads.txt.gz --ie1 ${jnc_path_prefix}_intron_exon_junction_counts.txt --ie2 ${jnc_path_prefix}_intron_exon_junction_counts.txt \
 -p ${outputDir}/${sample}/${sample}_${chr}/${sample}_${chr} --txt_db1 ${txt_db1} --txt_db2 ${txt_db2} --txt_db3 ${txt_db3} --jcn_seq_len ${jcn_seq_len} --sqlite_db_dir ${sqlite_db_dir}

