import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ground_truth")
    parser.add_argument("to_test")

    args = parser.parse_args()

    gt = set()
    with open(args.ground_truth) as fp:
        for line in fp:
            if line[0] != "@":
                fields = line.split()
                gt.add((fields[0], fields[3]))
    incorrect = 0
    total = 0
    unaligned = 0
    with open(args.to_test) as fp:
        for line in fp:
            if line[0] != "@":
                fields = line.split()
                if (int(fields[1]) & 4):
                    unaligned += 1
                elif (fields[0], fields[3]) not in gt:
                    incorrect += 1
                total += 1
    print(unaligned, incorrect, total)
