'''
sbcExp (c) University of Manchester 2017

sbcExp is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell
@description: SYNBIOCHEM experimental design factor extraction. 
'''
import re
import csv
import argparse
import openpyxl
import string
import logging
import datetime
import os
import pandas as pd


def timeStamp(arg):
    """ Create a unique time stamp for each run and log configuration"""
    ts = '{:%Y-%m-%d-%H-%M-%S}'.format(datetime.datetime.now())
    logfolder = os.path.join(arg.outFolder, 'log')
    if not os.path.exists(logfolder):
        os.mkdir(logfolder)
    logfile = os.path.join(logfolder, ts+'.log')
    logging.basicConfig(filename=logfile, level=logging.DEBUG)
    logging.info(ts)
    logging.info(arg.expFile)
    logging.info(arg.trackingSheet)
    logging.info(arg.doEInfo)
    logging.info(arg.levels)
    outfile = os.path.join(arg.outFolder,ts+'.csv')
    logging.info(outfile)
    return outfile

def arguments():
    parser = argparse.ArgumentParser('Extract DoE factors. Pablo Carbonell, SYNBIOCHEM, 2018')
    parser.add_argument('-doeFolder', default='/mnt/syno/shared/Designs/SBCDE00028/Design',
                        help='Experimental design folder')
    parser.add_argument('-doeFile', default=None,
                        help='Experimental design file (default  Design_id.j0)')
    parser.add_argument('-doeInfo', default=None,
                        help='Experimental design info file (default Design_id.ji0)')
    parser.add_argument('-levels', default=None,
                        help='Levels information for quantitative factors (id, description, level). Default: categorical factors')
    parser.add_argument('-outFolder', default=None,
                        help='Output folder (default: doeFolder)')
    return parser


def readLibrary(infoFile):
    dlib = {}
    for l in open(infoFile):
        m = l.rstrip().split('\t')
        sbcid = m[0]
        try:
            sbcid = str(int(re.sub('SBC', '', sbcid)))
        except:
            pass
        mm = []
        for x in m:
            x = re.sub('plasmid', 'l', x)
            x = re.sub('promoter', 'p', x)
            x = re.sub('gene', 'g', x)
            mm.append(x)
        dlib[sbcid] = mm[1:]
    return dlib

def contiguousGene(dlib, labels):
    """ Flag where 2 genes are contiguous """
    cont = {}
    for design in dlib:
        gprev = None
        for x in dlib[design]:
            m = x.split('_')
            if x.startswith('g'):
                if labels[x] == 'None':
                    # This is a special case of "empty gene"
                    # Keep gid as previous gene and continue
                    gid = gprev
                    continue
                gid = m[0]
                if gid in labels:
                    gid = labels[gid]
                if gprev != None:                    
                    key = gprev+'_'+gid
                    if key not in cont:
                        cont[key] = set()
                    cont[key].add(design)
                gprev = gid
            elif x.startswith('p') and labels[x] == 'None':
                gprev = gid
                continue
            else:
                gprev = None
    return cont
                
def promoterGene(dlib, labels):
    pg = {}
    for design in dlib:
        prom = None
        for x in dlib[design]:
            m = x.split('_')
            if x.startswith('p'):
                prom = labels[x]
            elif x.startswith('g') and prom is not None:
                gid = m[0]
                if gid in labels:
                    gid = labels[gid]
                key = 'p' +'_'+ re.sub('g', 'g', gid)
                if key not in pg:
                    pg[key] = {}
                pg[key][design] = prom
                prom = None
            else:
                prom = None
    return pg

def plasmids(dlib, labels):
    """ Extract plasmid type """
    pl = {}
    for design in dlib:
        plasm = dlib[design][0]
        pl[design] = labels[ plasm ]
    return pl

def genes(dlib, labels):
    """ Extract genes at each position """
    gl = {}
    for design in dlib:
        for x in dlib[design]:
            if x.startswith('g'):
                m = x.split('_')
                gid = 'g_'+m[0]
                if gid not in gl:
                    gl[gid] = {}
                gl[gid][design] = labels[ x ]
    return gl
                    
def parts(f):
    """ To do: Potentially we could reaplace the parts by some values,
    for instance for promoters. """
    px = {}
    lab = {}
    for l in open(f):
        m = l.split()
        if len(m) < 3:
            px[m[0]] = 0
        else:
            px[m[0]] = m[2]
        lab[m[0]] = m[1]
    return px, lab

def mapLibrary(dlib, equiv):
    """ Get the ids of each part. """
    eqlib = {}
    for line in open(equiv):
        construct = []
        line = line.rstrip()
        ids = line.split()
        plasmid = ids[0]
        l = 16
        ll = len(ids[1])
        for i in range(0, len(line)):
            if line[i:(i+ll)] == ids[1]:
                break
        while i < len(line):
            construct.append( line[i:(i+l)].rstrip() )
            i += l
        for j in range(0, len(dlib[plasmid])):
            try:
                if construct[j] == '':
                    construct[j] = 'None'
                if dlib[plasmid][j] in eqlib and construct[j] == 'None':
                    # do not update if it is not empty (see comment below)
                    continue
                eqlib[dlib[plasmid][j]] = construct[j]
            except:
                # This can happen because some parts are removed
                # in one of the constructs, prioritize assigning 
                # a value if somewhere happens in the library
                if dlib[plasmid][j] not in eqlib:
                    eqlib[dlib[plasmid][j]] = 'None'
            
    return eqlib


def dataExtraction(args):
    dlib = readLibrary(args.doeFile)
    dlib2 = mapLibrary(dlib, args.doeInfo)
    if args.levels is not None:
        px, lab = parts(args.levels)
    else:
        lab = dlib2
    cont = contiguousGene(dlib, lab)
    pg = promoterGene(dlib, lab)
    pl = plasmids(dlib, lab)
    gl = genes(dlib, lab)
    head = ['Design']
    head.extend(['pl']+sorted(gl)+sorted(pg)+sorted(cont))
    cw = []
    for design in sorted(dlib):
        row = [design]
        if design in pl:
            row.append(pl[design])
        else:
            row.append('None')
        for gen in sorted(gl):
            if design in gl[gen]:
                row.append(gl[gen][design])
            else:
                row.append('None')

        for prgen in sorted(pg):
            if design in pg[prgen]:
                row.append(pg[prgen][design])
            else:
                row.append('None')
        for cg in sorted(cont):
            if design in cont[cg]:
                row.append(1)
            else:
                row.append(0)
        cw.append( row )

    df = pd.DataFrame( cw, columns=head )
    return df

def mapFolder(doeFolder, doeFile=None, doeInfo=None, outFolder=None):
    parser = arguments()
    args = parser.parse_args(['-doeFolder', doeFolder])
    args.doeFolder = doeFolder
    if doeFile is not None:
        args.doeFile = doeFile
    if doeInfo is not None:
        args.doeInfo = doeInfo
    if outFolder is not None:
        args.outFolder = outFolder
    doMap(args)
    
def doMap(args):
    desname = os.path.basename(os.path.dirname(args.doeFolder))
    if args.doeFile is None:
        args.doeFile = os.path.join(args.doeFolder, desname+'.j0')
    if args.doeInfo is None:
        args.doeInfo = os.path.join(args.doeFolder, desname+'.ji0')
    if args.outFolder is None:
        args.outFolder = args.doeFolder
    outfile = os.path.join(args.doeFolder, desname+'_factors.csv')
    df = dataExtraction(args)
    df.to_csv( outfile, index= False )


if __name__ == '__main__':
    parser = arguments()
    args = parser.parse_args()
    doMap(args)
