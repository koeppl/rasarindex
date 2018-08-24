import sys
import argparse
from textwrap import fill

def read_first_seq(fasta_file):
    # returns (record_name, sequence)
    seq = ""
    record_name = ""
    with open(fasta_file) as fp:
        record_name = fp.readline().strip().split(">")[1]
        for line in fp:
            if line[0] == ">":
                break
            seq += line.strip()
    return (record_name, seq)

def chimeric_seq(seq1, seq2, segment_length):
    length = min((len(seq1), len(seq2)))
    new_seq = ""
    for x, i in enumerate(range(0, length, segment_length)):
        if x % 2 == 0:
            new_seq += seq1[i:i+segment_length]
        else:
            new_seq += seq2[i:i+segment_length]
    return new_seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1")
    parser.add_argument("fasta2")
    parser.add_argument("-s", "--segment_length",type=int, default=120000, help="""length of each partition of
            the resulting chimeric sequence. If segment_length >
            len(end_result)/2, then the first,larger segment will come from
            fasta1""")
    args = parser.parse_args()
    # by default, only extract the first record from each fasta file

    ref1, seq1 = read_first_seq(args.fasta1)
    ref2, seq2 = read_first_seq(args.fasta2)
    new_ref = "{}+{}".format(ref1, ref2)
    new_seq = chimeric_seq(seq1, seq2, args.segment_length)
    sys.stdout.write(">{}\n{}\n".format(new_ref, fill(new_seq, width=60)))
