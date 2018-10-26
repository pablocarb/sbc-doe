'''
brOligos (c) University of Manchester 2015

brOligos is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: Pablo Carbonell, SYNBIOCHEM
@description: Generates the list of bridging oligos for LCR assembly. Integrated within the general SYNBIOCHEM DoE pipeline
'''
 
from __future__ import print_function
import re
import argparse

def outInfo(message, outHandler=None, printInfo=True):
    try:
        outHandler.write(message+'\n')
    except:
        pass
    if printInfo:
        print(message)

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
    
def bridges(ct, ci, logFile=None):
    pool = set()
    ori = set()
    res = set()
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
               res.add(end) 
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
    outh = None
    if logFile is not None:
        try:
            outh = open(logFile, 'w')
        except:
            pass
    oriRes = False
    geneProm = set()
    geneOri = set()
    pool2 = set()
    outInfo('Broligos generator', outh)
    outInfo('(Ori, Res, Prom1, iProm, Gene) =  (%d, %d, %d, %d, %d)' % (len(ori), len(res),  len(pro1), len(prom), len(gene)), outh)
    n = 0
    # Pairing: Origin - Resistance
    for o in ori:
        # For resistance
        for r in res:
            # Take only one
            if not oriRes:
                pool2.add( (o,r) )
                oriRes = True
        
    outInfo('%d pairings: Origin - Resistance' % (len(pool2) - n,), outh)
    n += len(pool2)
    
    # Pairing: Resistance - Promoter
    for r in res:
        for p in pro1:
            pool2.add( (r,p) )
    outInfo('%d pairings: Resistance - First promoter' % (len(pool2) - n,), outh)
    n = len(pool2)

    for p in pro1:
        for g in gene:
            if g not in geneProm:
                pool2.add( (p,g) )
                geneProm.add( g )
    outInfo('%d pairings: First promoter - Gene' % (len(pool2) - n,), outh)
    n = len(pool2)
    
    for p in prom:
        for g in gene:
            pool2.add( (p,g) )
    outInfo('%d pairing: Intergenic promoters - Gene' % (len(pool2) - n,), outh)
    n = len(pool2)


    for g in gene:
        for p in prom:
            pool2.add( (g,p) )
    outInfo('%d pairings: Gene - Intergenic promoters' % (len(pool2) - n,), outh)
    n = len(pool2)

    # Pairing: Gene - Gene
    for g1 in gene:
        for g2 in gene:
            if g1 == g2:
                continue
            pool2.add( (g1, g2) )
            pool2.add( (g2, g1) )
    outInfo('%d pairings: Gene - Gene' % (len(pool2) - n,), outh)
    n = len(pool2)

    # Pairing: Gene - Origin
    for g in gene:
        for o in ori:
            if g not in geneOri:
                pool2.add( (g,o) )
                geneOri.add( g )
    outInfo('%d pairings: Gene - Origin' % (len(pool2) - n,), outh)
    n = len(pool2)

    outInfo('%d total pairings' % (n,), outh)

    try:
        outh.close()
    except:
        pass

    return pool, pool2

def arguments():
    parser = argparse.ArgumentParser(description='Bridge count, SYNBIOCHEM, 2018')
    parser.add_argument('doe', 
                        help='DoE file with ICE numbers')
    parser.add_argument('info', 
                        help='Info about factors')
    parser.add_argument('-outFile', 
                        help='Output file')
    parser.add_argument('-logFile', 
                        help='Log file')
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
    pool, pool2 = bridges(ct, ci, arg.logFile)
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
