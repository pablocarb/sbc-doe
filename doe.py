import numpy as np
from os import path
import pyRserve
import sys
import argparse

def construct(f):
    ct = []
    for l in open(f):
        m = l.rstrip().split('\t')
        if m[0].startswith('#'):
            continue
        factor = m[0]
        nlevels = int(m[1])
        try:
            positional = int(m[2])
        except:
            positional = 0
        ll = np.log2(nlevels)
        if ll != int(ll):
            raise 
        ct.append((factor, nlevels, positional))
    return ct

def getfactors(ct, permut=False):
    factors =[]
    nlevels = []
    npos = []
    i = 0
    for x in ct:
        f = x[0]
        l = x[1]
        try:
            p = x[2]
        except:
            p = 0
        if l>1:
            factors.append(f)
            nlevels.append(l)
        if p != 0:
            npos.append(x[0])
        i += 1
    return factors, nlevels, npos

# Assume promoters "p"  and count number of segments
def segments(libr):
    col = set()
    prom = set()
    for l in libr:
        x = ''
        pr = []
        for ll in l:
            m = ll.split('_')
            if ll.startswith('p'):
                if m[1] == '1':
                    if x != '':
                        pr.append(x)
                        col.add(x)
                    x = ''
            elif ll.startswith('g'):
                x += ll
        col.add(x)
        pr.append(x)
        prom.add('.'.join(sorted(pr)))

    if x != '':
        col.add(x)

    sta= {}
    for i in col:
        m = i.split('_')
        l = len(m)
        if l not in sta:
            sta[l] = 0
        sta[l] += 1
    return prom



parser = argparse.ArgumentParser(description='SBC-DeO. Pablo Carbonell, SYNBIOCHEM, 2016')
parser.add_argument('-p', action='store_true',
                    help='Full positional permutation (default: random latin square)')
parser.add_argument('-f', 
                    help='Input file with specifications')
parser.add_argument('-i', action='store_true',
                    help='Ignore segment calculations based on promoters')
arg = parser.parse_args()
if 'f' not in vars(arg):
    parser.print_help()
    sys.exit()
f = vars(arg)['f']
p = vars(arg)['p']
if not path.exists(f):
    parser.print_help()
    sys.exit()
wd = path.dirname(path.realpath(__file__))
conn = pyRserve.connect()
conn.r.source(path.join(wd, 'mydeo.r'))
ct = construct(f)
factors, nlevels, npos = getfactors(ct)
lat = None
if len(npos) > 0:
    factors.append('pos')
    if not p:
        lat = conn.r.permut(len(npos), ptype='latin')
    else:
        lat = conn.r.permut(len(npos), ptype='full')
    nlevels.append(len(lat))
if not p:
    designid = f
else:
    designid = f+'.full'
finfow = open(designid+'.info', 'w')
if not p:
    dinfo =  "SBC-DoE; Factors: %d; Levels: %d; Positional: %d [Latin square]" % (len(factors), np.prod(nlevels), len(npos))
else:
    dinfo = "SBC-DoE; Factors: %d; Levels: %d; Positional: %d [Full permutations]" % (len(factors), np.prod(nlevels), len(npos))
print(dinfo)
finfow.write(dinfo+'\n')
doe = conn.r.doe1(factors=np.array(factors), nlevels=np.array(nlevels), timeout=30)
librl = []
for des in range(0, len(doe)):
    ndes = {}
    n = 0
    # Read the design for each factor
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
    # Add positional factor
    if 'pos' in doe[des]['design'].keys:
        ndes['pos'] = np.array(doe[des]['design']['pos'])
    # Store designs
    fname = designid+'.d'+str(des)
    of = open(fname, 'w')
    libr = []
    for x in range(0, n):
        ll = []
        for y in ct:
            # Randomize permutations using a latin square
            fa = y[0]
            if fa in npos:
                perm = ndes['pos'][x]
                fa = npos[lat[perm-1][npos.index(fa)]-1]
            of.write("%s_%d\t" % (fa, ndes[fa][x],))
            ll.append("%s_%d" %  (fa, ndes[fa][x],))
        of.write('\n')
        libr.append(ll)
    of.close()
    librl.append(libr)
    if vars(arg)['i']:
        dinfor =  " Design %d; Model S^%d; Library size: %d" % (des, des+1, len(libr))
    else:
        dinfor = " Design %d; Model S^%d; Library size: %d; Segments: %d" % (des, des+1, len(libr), len(segments(libr)))
    print(dinfor)
    finfow.write(dinfor+'\n')
finfow.close()


# col = set()
# prom = set()
# for l in libr:
#     x = ''
#     pr = []
#     for ll in l:
#         m = ll.split('_')
#         if ll.startswith('p'):
#             if m[1] == '1':
#                 if x != '':
#                     pr.append(x)
#                     col.add(x)
#                 x = ''
#         elif ll.startswith('g'):
#             x += ll
#     col.add(x)
#     pr.append(x)
#     prom.add('.'.join(sorted(pr)))

# if x != '':
#     col.add(x)

# sta= {}
# for i in col:
#     m = i.split('_')
#     l = len(m) - 1
#     if l not in sta:
#         sta[l] = 0
#     sta[l] += 1
