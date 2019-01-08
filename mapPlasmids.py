#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 10:43:04 2019

@author: pablo
@content: Map plasmids ICE id using the csv list of entries.
@example: python mapPlasmids.py entries.csv SBCDE00040 SBCDE00040_export_mapping.csv
"""
import argparse
import pandas as pd
from datetime import datetime

def arguments():
    parser = argparse.ArgumentParser(description='Map plasmids. Pablo Carbonell, SYNBIOCHEM, 2019')
    parser.add_argument('entries', 
                        help='Entries csv file')
    parser.add_argument('design',  
                        help='Design ID')
    parser.add_argument('output',  
                        help='Ouput csv file')
    return parser

if __name__ == '__main__':
    parser = arguments()
    arg = parser.parse_args()
    df = pd.read_csv( arg.entries )
    plist = {}
    for ix in df.index:
        ptype = df.loc[ix,'type']
        if ptype != 'PLASMID':
            continue
        name = df.loc[ix,'name']
        try:
            sbcid, plid = name.split('_')
        except:
            continue
        part = df.loc[ix,'partId']
        ts = df.loc[ix,'creationTime']
        ts = datetime.fromtimestamp(float(ts)/1000)
        if sbcid == arg.design and plid.startswith('PL'):
            if name not in plist:
              plist[name] = (part, ts)
            else:
                t1 = plist[name][1]
                if ts > t1:
                    plist[name] = (part, ts)
    rows = []
    for p in sorted(plist):
        rows.append( (p, plist[p][0]) )
    out = pd.DataFrame(rows, columns=['Name','ICE'])
    out.to_csv(arg.output,index=False)