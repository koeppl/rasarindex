import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=argparse.FileType("r"), default=sys.stdin)
    args = parser.parse_args()

    name = ""
    seq = ""
    for line in args.fasta:
        if line[0] == ">":
            if seq != "":
                # flush
                try:
                    sys.stdout.write(name.strip() + "\n")
                    sys.stdout.write(seq.strip() + "\n")
                    sys.stdout.write("+\n")
                    sys.stdout.write("~"*len(seq) + '\n')
                except BrokenPipeError:
                    sys.stdout.close()
                    args.fasta.close()
                    exit(0)
                # set up new
                seq = ""
            name = "@" + line[1:]
        else:
            seq += line.strip()
    try:
        sys.stdout.write(name.strip() + "\n")
        sys.stdout.write(seq.strip() + "\n")
        sys.stdout.write("+\n")
        sys.stdout.write("~"*len(seq) + '\n')
    except BrokenPipeError:
        sys.stdout.close()
        args.fasta.close()
        exit(0)
