'''
Unique (c) University of Manchester 2015

Unique is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell
@description: SYNBIOCHEM design of experiments
'''
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
        try:
            pool = int(m[3])
        except:
            pool = 0
        # Degeneration: higher levels go to 1
        try:
            deg = int(m[4])
        except:
            deg = 0
#        ll = np.log2(nlevels)
#        if ll != int(ll):
#            raise 
        ct.append((factor, nlevels, positional, pool, deg))
    return ct


def convert_construct(xct):
    rid = {}
    ct = []
    for p in sorted(xct):
        levels = xct[p]['levels']
        comp = xct[p]['component']
        # This is an exception for promoters, need to be improved
        deg = 0
        if comp  == 'promoter':
            deg = len(levels) - len(filter( lambda h: h is None, levels))
            if len(levels) > deg:
                deg += 1
        elif comp == 'gene':
            deg = len(levels) - len(filter( lambda h: h is None, levels))
            if len(levels) ==  deg:
                deg = 0
        pos = xct[p]['positional']
        if pos is None:
            pos = 0
        cid = str(xct[p]['component'])+str(p)
        i = 1
        if comp == 'promoter' and len(levels) > len(filter( lambda h: h is None, levels)):
            rid[cid+'_'+str(i)] = None
            i += 1
        for x in range(0, len(levels)):
            if levels[x] is not None:
                rid[cid+'_'+str(i)] = levels[x]
                i += 1
        ct.append((cid, len(levels), str(pos), 0, deg))
        
    return ct, rid

def read_excel(e):
    import openpyxl
    wb = openpyxl.load_workbook(e)
    xl = wb.get_sheet_by_name(wb.get_sheet_names()[0])
    nl = xl.get_highest_row()
    fact = {}
    for r in range(1, nl):
        factor = int(xl.cell(row=r, column= 0).value)
        positional = xl.cell(row=r, column= 1).value
        component = xl.cell(row=r, column= 2).value
        part = xl.cell(row=r, column= 3).value
        if factor not in fact:
            fact[factor] = {'positional': positional,
                            'component': component,
                            'levels': []}
        if part == 'blank':
            part = None
        fact[factor]['levels'].append(part)
        
    return fact


def getfactors(ct, permut=False):
    factors =[]
    nlevels = []
    npos = []
    i = 0
    for x in ct:
        f = x[0]
        l = x[1]
        try:
            p = int(x[2])
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
def segments(libr, ct):
    col = set()
    prom = set()
    for l in libr:
        x = ''
        pr = []
        for i in range(0, len(l)):
            ll = l[i]
            cti = ct[i][1]
            m = ll.split('_')
            if ll.startswith('p'):
                newsegment = False
                # If promoter has only 1 level, we assume always on
                if cti == 1:
                    newsegment = True
                    plevel = 2 # level 1 means off, level 2 on
                # else, level 1 means promoter off
                else:
                    # Promoter 1 is always on (add one to the level, to be improved in a more general way)
                    if m[0] == 'p1':
                        newsegment = True
                        plevel = int(m[1]) + 1
                    elif m[1] != '1':
                        newsegment = True
                        plevel = int(m[1])
                if newsegment:
                    if x != '':
                        pr.append(x)
                        col.add(x)
#                    x = ''
                    x = str(plevel) # Put only the promoter level?
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
    return col
#    return prom


    

def save_design(design, ct, fname, lat, rid = None):
    ndes = {}
    n = 0
    # Read the design for each factor
    for x in ct:
        fact = x[0]
        if fact in design['design'].keys:
            ndes[fact] = np.array(design['design'][fact])
            n = len(ndes[fact])
        else:
            ndes[fact] = np.array([])
    for x in ct:
        fact = x[0]
        nlevels = x[1]
        dege = x[4]
        if dege > 0:
            for i in range(0, len(ndes[fact])):
                if ndes[fact][i] > dege:
                    ndes[fact][i] = 1
    for x in ct:
        fact = x[0]
        if len(ndes[fact]) == 0:
            ndes[fact] = np.repeat(1, n)
    # Add positional factor
    if 'pos' in design['design'].keys:
        ndes['pos'] = np.array(design['design']['pos'])
    # Store designs
    of = open(fname, 'w')
    libr = []
    libscr = []
    for x in range(0, n):
        ll = []
        screen = 1
        for y in ct:
            fa = y[0]
            le = y[1]
            po = y[2]
            pl = y[3]
            de = ndes[fa][x]
            # Screening pool?
            if pl > 0:
                if de > 1 or (le  == 1 and de > 0):
                    screen *= pl                        
            # Randomize permutations using a latin square
            faid = "%s_%d" % (fa, ndes[fa][x],)
            if fa in npos:
                perm = ndes['pos'][x]
                fa = npos[lat[perm-1][npos.index(fa)]-1]
                faid = "%s_%d" % (fa, ndes[fa][x],)
            if rid is not None:
                faid = rid[faid]
            if faid is not None:
                of.write("%s\t" % (faid,))
                ll.append("%s" %  (faid,))
        of.write('\n')
        libr.append(ll)
        if screen > 1:
            screen *= 3 # if screening a pool, multiply by 3
        libscr.append(screen)
    of.close()
    return libr, libscr

parser = argparse.ArgumentParser(description='SBC-DeO. Pablo Carbonell, SYNBIOCHEM, 2016')
parser.add_argument('-p', action='store_true',
                    help='Full positional permutation (default: random latin square)')
parser.add_argument('-f', 
                    help='Input file with specifications (old format)')
parser.add_argument('-e', 
                    help='Input file with specifications (excel format)')
parser.add_argument('-i', action='store_true',
                    help='Ignore segment calculations based on promoters')
parser.add_argument('-r', action='store_false',
                    help='No regular fractional factorial design')
parser.add_argument('-o', action='store_false',
                    help='No orthogonal array design')
parser.add_argument('-x', action='store_true',
                    help='Randomize orthogonal array design')
arg = parser.parse_args()
if 'f' not in vars(arg) and 'e' not in vars(arg):
    parser.print_help()
    sys.exit()
f = vars(arg)['f']
e = vars(arg)['e']
if f is not None:
    ct = construct(f)
    rid = None
else:
    if e is None or not path.exists(e):
        parser.print_help()
        sys.exit()
    xct = read_excel(e)
    ct, rid = convert_construct(xct)
p = vars(arg)['p']
wd = path.dirname(path.realpath(__file__))
conn = pyRserve.connect()
conn.r.source(path.join(wd, 'mydeo.r'))
factors, nlevels, npos = getfactors(ct)
lat = None
if len(npos) > 0:
    factors.append('pos')
    if not p:
        lat = conn.r.permut(len(npos), ptype='latin')
    else:
        lat = conn.r.permut(len(npos), ptype='full')
    nlevels.append(len(lat))
of = f
if of is None:
    of = e
if not p:
    designid = of
else:
    designid = of+'.full'
finfow = open(designid+'.info', 'w')
if not p:
    dinfo =  "SBC-DoE; Factors: %d; Levels: %d; Positional: %d [Latin square]" % (len(factors), np.prod(nlevels), len(npos))
else:
    dinfo = "SBC-DoE; Factors: %d; Levels: %d; Positional: %d [Full permutations]" % (len(factors), np.prod(nlevels), len(npos))
print(dinfo)
finfow.write(dinfo+'\n')
if vars(arg)['r']:
    doe1 = conn.r.doe1(factors=np.array(factors), nlevels=np.array(nlevels), timeout=30)
    for des in range(0, len(doe1)):
        fname = designid+'.d'+str(des)
        libr, libscr = save_design(doe1[des], ct, fname, lat, rid)
        if rid is not None:
            fname = designid+'.di'+str(des)
            libr, libscr = save_design(doe1[des], ct, fname, lat, rid=None)
        if vars(arg)['i']:
            dinfor =  " Design %d; Model S^%d; Library size: %d" % (des, des+1, len(libr))
        else:
            dinfor = " Design %d; Model S^%d; Library size: %d; Segments: %d; Screening size: %d" % (des, des+1, len(libr), len(segments(libr, ct)), np.sum(libscr))
        print(dinfor)
        finfow.write(dinfor+'\n')
if vars(arg)['o']:
    if vars(arg)['x']:
        seed = np.random.randint(1e6)
    else:
        seed = 100
    doe2 = conn.r.doe2(factors=np.array(factors), nlevels=np.array(nlevels),
                       timeout=30, seed=seed)
    for des in range(0, len(doe2)):
        fname = designid+'.oad'+str(des)
        libr, libscr = save_design(doe2[des], ct, fname, lat, rid)
        if rid is not None:
            fname = designid+'.oadi'+str(des)
            libr, libscr = save_design(doe2[des], ct, fname, lat, rid=None)
        if vars(arg)['i']:
            dinfor = " Orthogonal Array Design; Library size: %d" % (len(librs),)
        else:
            dinfor = " Orthogonal Array Design; Library size: %d; Segments: %d; Screening size: %d" % (len(libr), len(segments(libr, ct)), np.sum(libscr))
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
