#!/usr/bin/python3
import sys

def hamm(s1, s2):
    return [a != b for (a, b) in zip(s1, s2)].count(True)

def neighbours (pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A','C','G','T']
    neighbourhood = set()
    suffix_neighbours = neighbours(pattern[1:], d)
    for suffix_neighbour in suffix_neighbours:
        if hamm(pattern[1:], suffix_neighbour) < d:
            for base in 'ACGT':
                neighbourhood.add(base + suffix_neighbour)
        else:
            neighbourhood.add(pattern[0] + suffix_neighbour)
    return neighbourhood

i = 0
for i, x in enumerate(neighbours("CATTT", 3)):
    sys.stdout.write("@{}\n{}\n+\n{}\n".format(i, x, "~"*len(x)))

