#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 09:41:25 2019

@author: pablo
@content: Get ICE entries from a collection and output to a csv file
"""
import os
from synbiochem.utils import ice_utils, net_utils
import pandas as pd
import json, argparse


def arguments():
    parser = argparse.ArgumentParser(description='Get ICE entries. Pablo Carbonell, SYNBIOCHEM, 2019')
    parser.add_argument('outcsv', 
                        help='Ouput csv file')
    parser.add_argument('-collection', default='SHARED', 
                        help='Collection (default SHARED)')
    parser.add_argument('-limit', default=100000, 
                        help='Limit (default 100000)')
    parser.add_argument('-chunk', default=100, 
                        help='Chunk (default 100000)')
    return parser

def response2df( response ):
    resp = json.loads( response )
    df = pd.DataFrame( resp['data'] )
    return df, resp['resultCount']

if __name__ == '__main__':
    parser = arguments()
    arg = parser.parse_args()
    iceurl = 'https://ice.synbiochem.co.uk'
    client = ice_utils.ICEClient(iceurl,os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
    url = iceurl+ '/rest/collections/'+arg.collection+'/entries?limit='+str( arg.chunk )
    response = net_utils.get( url, headers=client._ICEClient__headers )
    df, total = response2df( response )
    while df.shape[0] < total and df.shape[0] < arg.limit:
        print(df.shape[0])
        url = iceurl+ '/rest/collections/'+arg.collection+'/entries?limit='+str( arg.chunk )+'&offset='+str( df.shape[0] )
        response = net_utils.get( url, headers=client._ICEClient__headers )    
        df1, total = response2df( response )
        df = pd.concat( [df, df1], ignore_index=True, sort=True)
    if 'Name' not in df.columns:
        # Use both namings to avoid issues
        try:
            df['Name'] = df['name']
            df['Part ID'] = df['partId']
        except:
            pass
    df.to_csv( arg.outcsv, index=False )
