s1 = "CATTT"

def hamm(s1, s2):
    h = 0
    for x, y in zip(list(s1), list(s2)):
        if x != y:
            h += 1
    return h

fa = ""
with open("datasets/seq.fa") as fp:
    fp.readline()
    for line in fp:
        fa += line.strip()
window_size = 5
dist = 3
stuff = set()
x = 0
for i in range(0, len(fa)-window_size+1):
    s2 = fa[i:i+window_size]
    if hamm(s1, s2) < dist:
        stuff.add(s2)
        x += 1

print(len(stuff))
print(x)
