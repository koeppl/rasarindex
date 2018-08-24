import sys

def read_pop(fname):
    """ retuns {sample: {"pop": x, "super_pop": y, ...}, ...}
    """
    samples = {}
    fp = open(fname)
    field_names = fp.readline().strip().split()
    for line in fp:
        fields = line.strip().split()
        samples[fields[0]] = dict(zip(field_names[1:], fields[1:]))
    return samples

# def read_sam(fname, haplotype=None):
#     """ note: removes multimapped queries
#     returns: {pos: set(refs)}, references
#     """
#     sam = {}
#     all_refs = set()
#     cur_query = ""
#     refs = set()
#     good = True
#     for line in open(fname):
#         fields = line.strip().split("\t")
#         query = fields[0]
#         ref = fields[2]
#         pos = int(fields[3])
#         if cur_query != query:
#             if cur_query != "" and good:
#                 # TODO: note: we're indexing by POSITION, not the query name
#                 if pos not in sam:
#                     if (haplotype and haplotype in cur_query) or (not haplotype):
#                         sam[pos] = {'refs': refs, "queries": [cur_query]}
#                 else:
#                     if (haplotype and haplotype in cur_query) or (not haplotype):
#                         sam[pos]['refs'].update(refs)
#                         sam[pos]['queries'].append(cur_query)
#                 all_refs = all_refs.union(refs)
#             elif not good:
#                 sys.stderr.write("skipping {} due to multimapping\n".format(cur_query))
#             cur_query = query
#             good = True
#             refs = set()
#         if ref in refs:
#             good = False
#         elif good:
#             refs.add(ref)
#     if cur_query != "" and good:
#         # TODO: note: we're indexing by POSITION, not the query name
#         if pos not in sam:
#             if (haplotype and haplotype in cur_query) or (not haplotype):
#                 sam[pos] = {'refs': refs, "queries": [cur_query]}
#         else:
#             if (haplotype and haplotype in cur_query) or (not haplotype):
#                 sam[pos]['refs'].update(refs)
#                 sam[pos]['queries'].append(cur_query)
#         all_refs = all_refs.union(refs)
#     elif not good:
#         sys.stderr.write("skipping {} due to multimapping\n".format(cur_query))
#     return sam, all_refs



def read_sam(fname, haplotype=None):
    all_refs = set()
    refs = set()
    good = True
    prev_query = ""
    prev_pos = 0
    sam = {}
    ps = set()
    for line in open(fname):
        fields = line.strip().split("\t")
        query = fields[0]
        ref = fields[2]
        # TODO: fix this
        flag = int(fields[1])
        if (flag & 4):
            continue
        pos = int(fields[3])
        if query != prev_query:
            #flush
            if prev_query != "" and len(ps) == 1:
                all_refs.update(refs)
                p = list(ps)[0]
                if p not in sam: 
                    if (haplotype and haplotype in prev_query) or (not haplotype):
                        sam[p] = {'refs':refs}
                else:
                    if (haplotype and haplotype in prev_query) or (not haplotype):
                        sam[p]['refs'].update(refs)
            elif len(ps) != 1:
                sys.stderr.write("skipping {} due to multimapping ({} occurences)\n".format(prev_query, len(ps)))
            # reset
            refs = set()
            good = True
            ps = set()
        refs.add(ref)
        prev_query = query
        ps.add(pos)
    # do the last line
    if prev_query != "" and len(ps) == 1:
        all_refs.update(refs)
        p = list(ps)[0]
        if p not in sam: 
            if (haplotype and haplotype in prev_query) or (not haplotype):
                sam[p] = {'refs':refs}
        else:
            if (haplotype and haplotype in prev_query) or (not haplotype):
                sam[p]['refs'].update(refs)
    elif len(ps) != 1:
        sys.stderr.write("skipping {} due to multimapping ({} occurences)\n".format(prev_query, len(ps)))
    return sam, all_refs
