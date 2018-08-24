#!/usr/bin/python

import sys
import argparse

def get_seeds(seq, k, s):
    seeds = []
    for i in range(0, len(seq)-k+1, s+1):
        seeds.append(seq[i:i+k])
    return seeds


def print_to_fastq(seqs, fp):
    for i, s in enumerate(seqs):
        fp.write("@seq_{}\n{}\n+\n{}\n".format(i, s, "~"*len(s)))

def print_to_pizza(seqs, fp):
    """ assume that all the sequences are of the same length """
    number = len(seqs)
    length = len(seqs[0])
    fp.write("# number={} length={} file= forbidden=\\n\\t\n".format(number, length))
    fp.write("{}\n".format("\n".join(seqs)))

if __name__ == "__main__":
    """ usage: python3 read_to_seeds.py <infile> -k <seed length> -s <seed_spacing> -f <output_format>

    <infile> : fastq file containing a single sequence
    <seed_lengh> : integer representing seed length (default 8)
    <seed_spacing> : integer describing how far to space seeds apart (default 0)
    <output_format> : [fastq, pizza] (default fastq)
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("-k", "--seed_length", type=int, default=8)
    parser.add_argument("-s", "--seed_spacing", type=int, default=0)
    parser.add_argument("-f", "--output_format", choices=["fastq", "pizza"], default="fastq")
    parser.add_argument("-if", "--input_format", choices=["fasta", "fastq"], default="fastq")

    args = parser.parse_args()

    seeds = []

    if args.input_format == "fastq":
        with open(args.infile) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 1:
                    read = line.strip()
                    seeds += get_seeds(read, args.seed_length, args.seed_spacing)

    elif args.input_format == "fasta":
        with open(args.infile) as fp:
            seq = ""
            for line in fp:
                if ">" in line:
                    if seq != "":
                        # start new seq
                        seeds += get_seeds( seq, args.seed_length, args.seed_spacing)
                    seq = ""
                else:
                    seq += line.strip()
            if seq != "":
                seeds += get_seeds(seq, args.seed_length, args.seed_spacing)
    
    if args.output_format == "fastq":
        print_to_fastq(seeds, sys.stdout)
    elif args.output_format == "pizza":
        print_to_pizza(seeds, sys.stdout)

