#! /usr/bin/env python3


import argparse
from statistics import median

def parse_args():
    parser = argparse.ArgumentParser(description='filters and converts the format of juncBase output to make it usable for downstream analyses')
    parser.add_argument('-d', '--drimtable', required=True,
                        help='drimseq formatted table - make sure sample headers are formatted as samplename_condition_batch')
    parser.add_argument('-o', '--outfile', required=True,
                        help='output file name')
    parser.add_argument('--cond1', required=True,
                        help='condition 1 to be tested - should be in sample names')
    parser.add_argument('--cond2', required=True,
                        help='condition 2 to be tested - should be in sample names')
    parser.add_argument('--min_event_reads_per_sample', type=int, default=15)
    parser.add_argument('--min_inc_reads_per_sample', type=int, default=5)
    parser.add_argument('--min_delta_psi', type=int, default=10)

    args = parser.parse_args()
    return args


def parse_drimtable(file):
    eventtocounts = {}
    header = None
    for line in open(file):
        line = line.rstrip().split('\t')
        if line[0] != '':
            counts = [int(x) for x in line[3:]]
            id = line[1]
            if id not in eventtocounts: eventtocounts[id] = [counts]
            else: eventtocounts[id].append(counts)
        else:
            header = line
    return header, eventtocounts

def filter_event(args, conditions, inclusion, exclusion):
    totexp = [inclusion[x] + exclusion[x] for x in range(len(inclusion))]
    if sum([x >= args.min_event_reads_per_sample for x in totexp]) >= 6:
        if sum([x >= args.min_inc_reads_per_sample for x in inclusion]) >= 3 \
                or sum([x >= args.min_inc_reads_per_sample for x in exclusion]) >= 3:
            psi = [inclusion[x] / totexp[x] if totexp[x] > 0 else 0 for x in range(len(inclusion))]
            c1 = [psi[x] for x in range(len(psi)) if conditions[x] == args.cond1]
            c2 = [psi[x] for x in range(len(psi)) if conditions[x] == args.cond2]
            delta = median(c2) - median(c1)
            if abs(delta) * 100 >= args.min_delta_psi:
                return True
    return False


if __name__ == "__main__":
    args = parse_args()
    header, eventtocounts = parse_drimtable(args.drimtable)

    out = open(args.outfile, 'w')
    out.write('\t'.join(header) + '\n')
    conditions = [x.split('_')[-2] for x in header[3:]]

    c = 0
    for id in eventtocounts:
        inclusion, exclusion = eventtocounts[id]
        if filter_event(args, conditions, inclusion, exclusion):
            out.write('\t'.join([str(c), id, 'inclusion_' + id] + [str(x) for x in inclusion]) + '\n')
            c += 1
            out.write('\t'.join([str(c), id, 'exclusion_' + id] + [str(x) for x in exclusion]) + '\n')
            c += 1

