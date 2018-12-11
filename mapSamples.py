'''
mapSamples (c) University of Manchester 2018

mapSamples is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Map samples into their DoEs ids (and generate factors table)
'''
import os, re, argparse, csv
import pandas as pd
from synbiochem.utils import ice_utils

def arguments():
    parser = argparse.ArgumentParser(description='MapSamples. Pablo Carbonell, SYNBIOCHEM, 2018')
    parser.add_argument('dts', 
                        help='Data tracking sheet')
    parser.add_argument('-outFile',  default=None,
                        help='Output file (default: dtsFolder/samples.csv)')
    return parser

def samples(dts, sheetName='Samples'):
    try:
        df = pd.read_excel(dts, sheetName)
    except:
        df = None
    return df

def mapPlasmids(df, column='Strain ID', mapColumn='Plasmid Name'):
    ids = set()
    for x in df.loc[:,column]:
        try:
            iceid = int(x)
            ids.add( iceid )
        except:
            continue
    client = ice_utils.ICEClient("https://ice.synbiochem.co.uk",os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
    plmap = {}
    for ice in ids:
        try:
            entry = client.get_ice_entry( ice )
            name = entry.get_name()
            plmap[ ice ] = name
        except:
            continue
    for i in range(0, df.shape[0]):
        iceid = df.loc[i, column]
        if iceid in plmap:
            df.loc[i, mapColumn] = plmap[ iceid ]
    return df

def outputSamples(df, outfile):
    df.to_csv( outfile )

def readDesign( dfile, des={} ):
    with open(dfile) as h:
        for line in h:
            row = line.rstrip().split('\t')
            des[ row[0] ] = row[1:]

def addDesignColumns( df, designsFolder='/mnt/syno/shared/Designs', ext='.j0', mapColumn='Plasmid Name'):
    """ Create additional columns with the design combinations """
    """ TO DO: transform positional parameters into: promoters in front of genes, gene variants, and pairings  """
    designs = set()
    for plasmid in df[mapColumn].unique():
        try:
            if plasmid.startswith('SBCDE'):
                designs.add( plasmid.split('_')[0] )
        except:
            continue
    des = {}
    for desi in designs:
        dfile = os.path.join( designsFolder, desi, 'Design', desi+ext )
        readDesign( dfile, des )
    for i in range(0, df.shape[0]):
        plasmid = df.loc[i,mapColumn]
        if plasmid in des:
            for j in range(0, len( des[plasmid] ) ):
                df.loc[i,'pos'+str(j)] = des[plasmid][j]
    
    
    

if __name__ == '__main__':
    parser = arguments()
    arg = parser.parse_args()
    if arg.outFile is None:
        arg.outFile = os.path.join( os.path.dirname( arg.dts ), 'samples.csv' )
    df = samples( arg.dts )
    df = mapPlasmids( df )
    addDesignColumns( df )
    outputSamples( df, arg.outFile )
