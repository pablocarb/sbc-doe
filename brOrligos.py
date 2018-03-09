# Generate list of bridging oligos: to be integrated with the general DoE software
from __future__ import print_function
import re
import argparse

def readFiles(doe, info):
    ct = {}
    ci = {}
    with open(doe) as h:
        for line in h:
            row = line.rstrip().split('\t')
            cid = row[0]
            ct[cid] = row[1:]
    with open(info) as h:
        for line in h:
            row = line.rstrip().split('\t')
            cid = row[0]
            ci[cid] = []
            for part in row[1:]:
                if part.startswith('promoter') and re.match('.*_[3,4]$', part):
                        continue
                ci[cid].append( part )
    return ct, ci
    
def bridges(ct, ci):
    pool = set()
    ori = set()
    pro1 = set()
    prom = set()
    gene = set()
    for x in sorted(ct):
        start = None
        if len(ct[x]) != len(ci[x]):
            import pdb
            pdb.set_trace()
        for i in range(0, len(ct[x])):
            part = ci[x][i]
            ice = ct[x][i]
            end = ice
            if part.startswith('resistance'):
                continue
            if part.startswith('origin'):
                ori.add( end )
            elif part.startswith('promoter3'):
                pro1.add( end )
            elif part.startswith('promoter'):
                prom.add( end )
            elif part.startswith('gene'):
                gene.add( end )
            if start is not None:
                pool.add( tuple( sorted( (start, end) ) ) )
            start = end
    pool2 = set()
    for o in ori:
        for p in pro1:
            pool2.add( (o,p) )
    for p in pro1 | prom:
        for g in gene:
            pool2.add( (p,g) )
            if p not in pro1:
                pool2.add( (g,p) )
    for g1 in gene:
        for g2 in gene:
            if g1 == g2:
                continue
            pool2.add( (g1, g2) )
            pool2.add( (g2, g1) )
    return pool, pool2

def arguments():
    parser = argparse.ArgumentParser(description='Bridge count, SYNBIOCHEM, 2018')
    parser.add_argument('doe', 
                        help='DoE file with ICE numbers')
    parser.add_argument('info', 
                        help='Info about factors')
    parser.add_argument('-outFile', 
                        help='Output file')
    return parser


def command_line(parser, args=None):
    if args is None:
        arg = parser.parse_args()
    else:
        arg = parser.parse_args(args)        
    return arg

def run_bro(args=None):
    parser = arguments()
    arg = command_line(parser, args)
    ct, ci =readFiles(arg.doe, arg.info)
    pool, pool2 = bridges(ct,ci)
    if arg.outFile is None:
        for x in sorted(pool2):
            print( '\t'.join( (x[0], x[1], '_'.join(x)) ) )
    else:
        try:
            with open(arg.outFile, 'w') as h:
                for x in sorted(pool2):
                    h.write( '\t'.join( (x[0], x[1], '_'.join(x)) )+'\n' )
        except:
            print( 'Error' )
                
    

if __name__ == "__main__":
    run_bro()
