import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("f", default=sys.stdin, type=argparse.FileType("r"))
    args = parser.parse_args()
    current_field = ""
    alignments = ""
    refs = set()
    print_this = True
    for line in args.f:
        fs = line.strip().split('\t')
        if current_field == "":
            current_field = fs[0]
        if fs[0] != current_field:
            # flush
            if print_this:
                sys.stdout.write(alignments)
            # reset
            alignments = ""
            refs = set()
            print_this = True
            current_field = fs[0]
        if fs[1] in refs:
            print_this = False
        else:
            refs.add(fs[1])
            alignments += line
    if print_this:
        sys.stdout.write(alignments)
