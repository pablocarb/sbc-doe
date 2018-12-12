# A test in order to compress a library containing odd levels
# Bad idea: log is not regularly spaced
import sys

equiv = {'4': '3'}

f = sys.argv[1]
for l in open(f):
    m = l.rstrip().split('\t')
    mm = []
    for x in m:
        a, b = x.split('_')
        if a in ['p2', 'p3', 'p4']:
            if b in equiv:
                b = equiv[b]
        mm.append('_'.join([a, b]))
    print '\t'.join(mm)
            
