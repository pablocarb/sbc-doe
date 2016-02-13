# Graphical representation of the library
import sys
from os import system
from glob import glob
import re


f = sys.argv[1]
i = 0
gl = []
for l in open(f):
    m = l.rstrip().split('\t')
    for x in m:
        v = x.split('_')
        if v[0] not in gl:
            gl.append(v[0])

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
            if v[0] == 'p1' or v[1] != '1':
                if v[0] != 'p1':
                    ow.write('t t%s\n' % (v[0][1],))
                ow.write('p %s\n' % (v[0],))
            else:
                continue
        else:
            ow.write('c %s %d\n' % (v[0],gl.index(v[0])+1))
    ow.write('t tx\n')
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
