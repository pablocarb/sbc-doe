# A test in order to compress a library containing odd levels
# Bad idea: log is not regularly spaced
import numpy as np

for l in open('160203_designs/test3.txt.d0'):
    m = l.rstrip().split('\t')
    mm = []
    for x in m:
        a, b = x.split('_')
        if a in ['p2', 'p3', 'p4']:
            b = str(int(np.log2(int(b)-1)+1))
        mm.append('_'.join([a, b]))
    print '\t'.join(mm)
            
