# Graphical representation of the library using pigeon cad
import sys
from os import system
from glob import glob
import re

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
        try:
            pool = int(m[3])
        except:
            pool = 0
        # Degeneration: in principle ignored
        try:
            deg = int(m[4])
        except:
            deg = 0
        ct.append((factor, nlevels, positional, pool, deg))
    return ct



# Read library
f = sys.argv[1]
if len(sys.argv)>2:
    dinfo = sys.argv[2]
    ct = construct(dinfo)
i = 0
gl = []
gl1 = []
# Count how many levels each factor has
count = {}
for l in open(f):
    m = l.rstrip().split('\t')
    for x in m:
        v = x.split('_')
        if v[0] not in gl:
            gl.append(v[0])
        if v[0] not in count:
            count[v[0]] = set()
        count[v[0]].add(x)
        if x not in gl1:
            gl1.append(x)

fl = []
for l in open(f):
    m = l.rstrip().split('\t')
    fn = f+'.pcad'+str(i)
    fl.append(fn)
    ow = open(fn, 'w')
    m = l.split()
    for x in m:
        v = x.split('_')
        if x.startswith('p'):
            # p1: assumes promoter absence
            if v[0] == 'p1' or v[1] != '1':
                if v[0] != 'p1':
                    ow.write('t\n')
#                ow.write('p %s %d\n' % (v[0],gl.index(v[0])+1))
# For promoters, we just give promoter number and color it accordingly
                ow.write('p p%s %d\n' % (v[1],int(v[1])*2+2))
            else:
                continue
        else:
#            ow.write('c %s %d\n' % (v[0],gl.index(v[0])+1))
            if len(count[v[0]]) > 1:
                   ow.write('c %s %d\n' % (x,gl1.index(x)+1))
            else:
                   ow.write('c %s %d\n' % (v[0],gl1.index(x)+1))                
    ow.write('t\n')
    ow.write('# Arcs\n')
    ow.close()
    i += 1

ofl = []
for pcad in fl:
    of = pcad+'.png'
    ofl.append(of)
    cmd = 'perl piget.pl '+pcad+' '+of
    system(cmd)

cmd = 'convert '+' '.join(ofl)+' -append '+f+'.png'
system(cmd)
