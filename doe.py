import numpy as np
from os import path
import pyRserve

def construct(f):
    ct = []
    for l in open(f):
        m = l.rstrip().split('\t')
        factor = m[0]
        nlevels = int(m[1])
        ll = np.log2(nlevels)
        if ll != int(ll):
            raise 
        ct.append((factor, nlevels))
    return ct

def getfactors(ct):
    factors =[]
    nlevels = []
    for x in ct:
        f = x[0]
        l = x[1]
        if l>1:
            factors.append(f)
            nlevels.append(l)
    return factors, nlevels

f = 'construct.txt'
wd = path.dirname(path.realpath(__file__))
conn = pyRserve.connect()
conn.r.source(path.join(wd, 'mydeo.r'))
ct = construct(f)
factors, nlevels = getfactors(ct)
doe = conn.r.doe1(factors=np.array(factors), nlevels=np.array(nlevels), timeout=30)
for des in range(0, len(doe)):
    ndes = {}
    n = 0
    for x in ct:
        fact = x[0]
        if fact in doe[des]['design'].keys:
            ndes[fact] = np.array(doe[des]['design'][fact])
            n = len(ndes[fact])
        else:
            ndes[fact] = np.array([])
    for x in ct:
        fact = x[0]
        if len(ndes[fact]) == 0:
            ndes[fact] = np.repeat(1, n)
    of = open(f+'.d'+str(des), 'w')
    for x in ct:
        of.write(x[0]+'\t')
    of.write('\n')
    for x in range(0, n):
        for y in ct:
            of.write("%d\t" % (ndes[y[0]][x],))
        of.write('\n')
    of.close()
