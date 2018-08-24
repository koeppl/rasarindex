import sys
import argparse

        
def binary_search(tree_arr, val):
    right = int(len(tree_arr)/2)
    if right == 0:
        return tree_arr[0]
    left = right - 1
    if val >= tree_arr[left] and val < tree_arr[right]:
        return tree_arr[left]
    elif val < tree_arr[left]:
        return binary_search(tree_arr[:left+1], val)
    else: # elif val >= tree_arr[right]:
        return binary_search(tree_arr[right:], val)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("offsets")
    parser.add_argument("faidx")
    parser.add_argument("--sam", "-s", action="store_true", help="specify if input file is in SAM format")
    args = parser.parse_args()

    chrm_offsets = []
    offset_to_chrm = {}

    with open(args.faidx) as fp:
        # chr length offset
        cur_off = 0
        for line in fp:
            fields = line.split()
            chrm = fields[0]
            length = int(fields[1])
            chrm_offsets.append(cur_off)
            offset_to_chrm[cur_off] = chrm
            cur_off += length 

    if args.offsets == "-": 
         fp = sys.stdin
    else:
        fp = open(args.offsets)
    if args.sam:
        for line in fp:
            fields = line.split("\t")
            flag = int(fields[1])
            if not (flag & 4):
                offset = int(fields[3])
                chrm_offset = binary_search(chrm_offsets, offset)
                real_offset = offset - chrm_offset
                fields[3] = str(real_offset)
                fields[2] = offset_to_chrm[chrm_offset]
            else:
                fields[3] = "*"
                fields[2] = "*"
            try:
                sys.stdout.write('\t'.join(fields))
            except BrokenPipeError:
               sys.stdout.close()
               exit()
    else:
        aln_offsets = map(int, fp.readlines())
        # TODO: does r-index do 1-based or 0-based?
        for offset in aln_offsets:
            # figure out which chromosome the offset corresponds to
            chrm_offset = binary_search(chrm_offsets, offset)
            real_offset = offset - chrm_offset
            try:
                # sys.stdout.write("{} {}\n".format(offset, chrm_offset))
                sys.stdout.write("{}\t{}\n".format(offset_to_chrm[chrm_offset], real_offset))
            except BrokenPipeError:
                sys.stdout.close()
                exit()
