# Provide a good example about segments

f = 'segments/s1.txt.d0'
prl = set()
asl = []
n = 0
for l in open(f):
    m = l.rstrip().split('\t')
    pr = []
    assembly = ''
    for x in m:
        v = x.split('_')
        if x.startswith('p'):
           if v[0] == 'p1':
               assembly = 'p' + str(int(v[1])+1)
           elif v[1] != '1':
               pr.append(assembly)
               assembly = 'p'+v[1]
        else:
            assembly += v[0]

    pr.append(assembly)
    asl.append(pr)
    print ' '.join(pr)
    for p in pr:
        prl.add(p)
    if len(pr) > n:
        n = len(pr)
print len(asl), len(prl)
prll = list(prl)
ow = open(f+'.seg', 'w')
for p in asl:
    v = []
    for i in range(0, n):
        if len(p) > i:
            c = p[i]
            ng = len(c.split('g')) -1
            for j in range(0, ng):
                v.append(prll.index(c)+1)
    ow.write(' '.join(map(str, v))+'\n')
ow.close()
