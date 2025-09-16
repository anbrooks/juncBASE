Core JuncBase functions are described in the MANUAL. Recommended postprocessing is described below:

If I run:
 - python3 juncbase_filteroutput_0925.py -j DATASET_AS_exclusion_inclusion_counts_lenNorm.tsv -o DATASET --gtf ANNOT.gtf

It will:
 - Remove jcn_only and novel intron_retention events
 - Re-annotate all events with genes, increasing the amount of events with unambiguous genes
 - Remove events that don’t have non-zero counts in at least 3 samples
 - Remove all overlapping events of the same type, so that each locus is represented by only the most interesting event (determined by median counts * stdev(PSI))
 - Generate a compressed, unique key for each event based on the most relevant junction/exon coordinates
 - Output both a drimseq formatted table using that key and a PSI table with all of the original juncbase info + that key as an extra column, so it’s easy to refer back to the full splicing event info

I will get:
 - DATASET.drimformat.tsv
 - DATASET.jb.dedup.psi.tsv
 - DATASET.dedup.bed

```
usage: juncbase_filteroutput_0925.py [-h] -j JUNCBASE -o OUTPUT [--gtf GTF]
filters and converts the format of juncBase output to make it usable for downstream analyses
options:
 -h, --help      show this help message and exit
 -j JUNCBASE, --juncbase JUNCBASE
            juncbase table with full event info and inclusion;exclusion counts
 -o OUTPUT, --output OUTPUT
            output prefix
 --gtf GTF       [optional] annotated gtf for improving gene definitions
```


If you have a large number of samples and want to do statisitcal analysis based on the PSI, you can stop here. If you want to run DRIMseq on this data, use the next step.

The file in the previous step can contain any number of samples, but before passing it into the next step, I recommend filtering to just the two conditions you want to use to run drimseq. You may also need to rename the columns to the correct format (sample_condition_batch)

Right before running DRIMseq, you can run:
 - python3 juncbase_drim_prefilter.py -d DATASET.drimformat.mysamples.tsv --cond1 WT --cond2 MUT -o DATASET.drimformat.mysamples.filtered.tsv

This will:
 - Remove events with low counts or low delta PSI (you can tweak what the exact parameters are)

```
usage: juncbase_drim_prefilter.py [-h] -d DRIMTABLE -o OUTFILE --cond1 COND1 --cond2 COND2
                 [--min_event_reads_per_sample MIN_EVENT_READS_PER_SAMPLE]
                 [--min_inc_reads_per_sample MIN_INC_READS_PER_SAMPLE] [--min_delta_psi MIN_DELTA_PSI]
filters drimseq formatted tables by PSI and counts
options:
 -h, --help      show this help message and exit
 -d DRIMTABLE, --drimtable DRIMTABLE
            drimseq formatted table - make sure sample headers are formatted as samplename_condition_batch
 -o OUTFILE, --outfile OUTFILE
            output file name
 --cond1 COND1     condition 1 to be tested - should be in sample names
 --cond2 COND2     condition 2 to be tested - should be in sample names
 --min_event_reads_per_sample MIN_EVENT_READS_PER_SAMPLE
 --min_inc_reads_per_sample MIN_INC_READS_PER_SAMPLE
 --min_delta_psi MIN_DELTA_PSI
```

