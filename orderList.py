#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 14:12:04 2018

@author: pablo

Void names for some ice
"""
import os,re
from synbiochem.utils import ice_utils
import pandas as pd
from datetime import datetime
import numpy as np
import csv


f = '/mnt/SBC1/data/ice/SBCDE00022.csv'
#ftest = '/mnt/SBC1/data/Biomaterials/Designs/SBCDE00040/Design/SBCDE00040_export.csv'
#f = ftest
desid = 'SBCDE00022'
#desidtest = 'test'
#desid = desidtest
outFolder = os.path.join('/mnt/SBC1/data/Biomaterials/Designs',desid,'Design')
of1 = os.path.join( outFolder, desid+'_export.csv')
of2 = os.path.join( outFolder, desid+'_export_mapping.csv')
client = ice_utils.ICEClient(iceurl,os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
with open(of1,'w') as h1, open(of2, 'w') as h2:
    head1 = ['Part ID','Name','Description','Sequence','Type']
    head2 = ['Name','ICE']
    cw1 = csv.writer(h1, lineterminator='\n')
    cw2 = csv.writer(h2, lineterminator='\n')
    cw1.writerow( head1 )
    cw2.writerow( head2 )
    # Read the list of parts of the design
    df = pd.read_csv(f)
    parts = set()
    for i in df.index:
        part_id = df.loc[i,'Part ID']
        if part_id in parts:
            continue
        entry = client.get_ice_entry( part_id )
        pType = entry.get_metadata()['type']
        if pType == 'PLASMID':
            design = entry.get_parameters()['Design id']
            name = entry.get_parameters()['Full name']
            cw2.writerow( [name, part_id] )
        else:
            pType = entry.get_parameters()['Type']
            if pType != 'DOMINO':
                continue
            name = entry.get_metadata()['name']
            items = name.split(' ')
            for i in range(0,len(items)):
                if items[i].startswith('SBC'):
                    parts.add(items[i])
            description = entry.get_metadata()['shortDescription']
            seq = entry.get_seq()
            cw1.writerow( [part_id, name, description, seq, pType] )
    for part_id in sorted(parts):
        entry = client.get_ice_entry( part_id )
        pType = ''
        for x in entry.get_metadata()['parameters']:
            if x['name'] == 'Type':
                pType = x['value']
        name = entry.get_metadata()['name']
        description = entry.get_metadata()['shortDescription']
        seq = entry.get_seq()
        cw1.writerow( [part_id, name, description, seq, pType] )
        
