# Script to convert our DoE design into data that can be read in JMP
from __future__ import print_function
import csv
import numpy as np

de = '/mnt/syno/shared/Designs/SBCDE00001/SBCDE00001.oadi0'
cv = csv.writer(open('/mnt/syno/shared/Pablo/jmp/SBCDE00001.csv', 'w'))
cv.writerow(['plasmid', 'gene1', 'promoter1', 'gene2', 'promoter2', 'gene3', 'promoter3', 'gene4', 'titer'])
for replica in range(0,3):
    for l in open(de):
        dd = []
        m = l.rstrip().split('\t')
        strength = 0.0
        for x in m:
            n = x.split('_')
            if x.startswith('plasmid'):
                dd.append(n[1])
            elif x.startswith('gene'):
                dd.append(n[0])
            elif x.startswith('promoter'):
                ps = int(n[1])-1
                strength += ps
                dd.append(str(int(n[1])-1))
        dd.append(2*np.random.sample()+strength*(0.8+0.2*np.random.sample()))
        cv.writerow(dd)
