import argparse
import re
import sys

softclipped_re = re.compile(r"^(\d+)S")
cigar_match_re = re.compile(r"(\d+)M")

def cigar_to_offset(cigar):
    # check if there's an S after the first string of numbers
    is_se = re.match(softclipped_re, cigar)
    # case 1: softclipping from left
    if (is_se):
        return int(is_se.groups(1)[0])
    # case 2: no softclipping
    else:
        return 0

def cigar_to_size(cigar):
    matches = re.findall(cigar_match_re, cigar)
    size = 0;
    if matches:
        for m in matches:
            size += int(m)
    return size

def min_aln_length(fields, length):
    return (fields[5] != "*" and cigar_to_size(fields[5]) >= length)

def flag(fields, tag):
    return (int(fields[1]) & tag) != 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file", type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("--min_aln_length", type=int, help="minimum alignment length (calculated from CIGAR string")
    parser.add_argument("--flag", type=int, help="SAM flag")
    # TODO: add other filters
    args = parser.parse_args()

    for line in args.sam_file:
        if line[0] != "@":
            fields = line.strip().split("\t")
            # if no options specified, to_print should evaluate to True
            to_print = (
                    ((not args.flag) or flag(fields, args.flag)) and
                    ((not args.min_aln_length) or min_aln_length(fields, args.min_aln_length))
                    ) # TODO add other filtering statements, stringed by AND
            if to_print:
                sys.stdout.write(line)
        else:
            sys.stdout.write(line)
