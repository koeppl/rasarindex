#!/bin/env python3
import sys
from utils import read_sam, read_pop
from collections import Counter
import argparse

if __name__ == "__main__":
    # sam = read_sam("alignments/HG00142.6.mhc.100bp.0.1x.1.sam", ".1")
    parser = argparse.ArgumentParser()
    parser.add_argument("sam")
    parser.add_argument("-t", "--haplotype", default=None)
    parser.add_argument("-d", "--default_ref")
    parser.add_argument("--no_read_len", default=False, action='store_true',
    help=""""use this option when the query lengths are informed by the
    previous query position. Make sure you use seq_length properly with this
    option""")
    parser.add_argument("-l", "--read_len", type=int, default=100)
    parser.add_argument("-N", "--seq_len", type=int, default=6000001)
    args = parser.parse_args()

    # read the sam file, scan and create running intersections
    # returns {position: {reference:set(), query:list()}}
    sam, all_refs = read_sam(args.sam, args.haplotype)
    # TODO: filter sam by ALIGNMENT SIZE
    sorted_positions = sorted(sam)
    running_ref_set = all_refs.copy()
    running_position_list = list()
    clusters = []
    ref_counter = Counter()
    for pos in sorted_positions:
        if pos == -1: #TODO: is this okay?
            continue
        ref_set = set.intersection(running_ref_set, sam[pos]['refs'])
        if len(ref_set) == 0:
            # print(running_position_list, running_ref_set)
            clusters.append((running_position_list, running_ref_set))
            ref_counter.update(running_ref_set)
            ref_set = sam[pos]['refs'].copy()
            running_position_list = []
        running_ref_set = ref_set
        running_position_list.append(pos)
    clusters.append((running_position_list, running_ref_set))
    ref_counter.update(running_ref_set)
    
    # print out faidx statements here
    refs_in_order = list(zip(*list(ref_counter.most_common())))[0]
    if args.default_ref:
        default_ref = args.default_ref
    else:
        default_ref = refs_in_order[0]
    end = 0
    for i, (positions, references) in enumerate(clusters):
        if args.no_read_len:
            assert(positions[0] - end == 0)
        # choose the reference
        if default_ref in references:
            ref = default_ref
        else:
            for r in refs_in_order:
                if r in references:
                    ref = r
                    break
        if args.no_read_len:
            # extend end to the next cluster's beginning if no read length is specified
            if i < len(clusters)-1:
                cur_end = clusters[i+1][0][0] 
            else:
                cur_end = args.seq_len
        else: 
            # extend end to read_len+last_position
            cur_end = positions[-1]+args.read_len
        if positions[0] - end > 0:
            # fill in the previous gap
            faidx_string = "{}:{}-{}".format(default_ref, end+1, positions[0])
            sys.stdout.write(faidx_string+"\n")
            faidx_string = "{}:{}-{}".format(ref, positions[0]+1, cur_end)
            sys.stdout.write(faidx_string+"\n")
        else:
            faidx_string = "{}:{}-{}".format(ref, end+1, cur_end)
            sys.stdout.write(faidx_string+"\n")
        end = cur_end
    # last bit
    if not args.no_read_len:
        faidx_string = "{}:{}-{}".format(default_ref, end+1, args.seq_len)
        sys.stdout.write(faidx_string+"\n")
