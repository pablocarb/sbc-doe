"""
After running mapSamples:
%run mapSamples.py /mnt/SBC1/data/stress2/dts -outFolder /mnt/SBC1/data/stress2/learn -iceEntries ice20190816.cs

Run this script in order to get the predicted best plasmids for the next round. Two possibilities:
- Select mimimum full factorial for top plasmids (12 combinations)
- Select top plasmids (12 by default)

To do: accept gene rearrangement
"""
import argparse
import os
import re
import pandas as pd
import numpy as np


def arguments():
    parser = argparse.ArgumentParser(description='Predict top constructs. Pablo Carbonell, SYNBIOCHEM, 2019')
    parser.add_argument('plasmidFile', 
                        help='Plasmid file contains numbered columns indicating the top combinations and their role. Template:1.origin,2.resistance,3.promoter,4.gene, etc... It might be either automatically generated through the mapSample.py script (learn analysis) or manually created.')
    parser.add_argument('-n', default=None,
                        help='Number of top plasmids')
    parser.add_argument('-s', default=None,
                        help='Maximum number of combinations')
    parser.add_argument('-designsFolder', default='/mnt/syno/shared/Designs',
                        help='Folder that contains the design templates')
    return parser

if __name__ == '__main__':
    parser = arguments()    
    arg = parser.parse_args()
    df = pd.read_csv(arg.plasmidFile)
    if arg.s is not None:
        """ Redefine the DoE """
        for i in np.arange(df.shape[0]):
            pos = []
            for col in df.columns:
                try:
                    int( col.split('.')[0] )
                except:
                    break
                pos.append( df.loc[0:i,col].unique() )
            space = np.prod([len(x) for x in pos])
            print(i,space)
            if space > int(arg.s):
                break
        newdes = {}
        for col in df.columns:
            try:
                posi = int( col.split('.')[0] )
            except:
                break
            newdes[posi] = set(df.loc[0:i-1,col])

        des = os.path.basename(arg.plasmidFile).split('_')[0]
        doeFile = os.path.join(arg.designsFolder,
                               des,
                               'Design',
                               'DoE_'+re.sub('SBCDE', 'SBC', des)+'.xlsx')
        if os.path.exists(doeFile):
            doeinfo = pd.read_excel( doeFile )
            newdoe = doeinfo.copy()
            for i in doeinfo.index:
                row = doeinfo.iloc[i,:]
                pos = int(row['DoE position'])
                part = row['Part number']
                if pos in newdes and part in newdes[pos]:
                    continue
                else:
                    newdoe = newdoe.drop(i)
        newdoe.to_excel(os.path.join(os.path.dirname(arg.plasmidFile),
                        re.sub('.csv','',os.path.basename(arg.plasmidFile))+'_newDoE.xlsx'),
                        index=False)
    else:
        des = os.path.basename(arg.plasmidFile).split('_')[0]
        doeFile = os.path.join(arg.designsFolder,
                               des,
                               'Design',
                               'DoE_'+re.sub('SBCDE', 'SBC', des)+'.xlsx')
        
        doeinfo = pd.read_excel( doeFile )
        levels = {}
        for i in doeinfo.index:
            row = doeinfo.iloc[i,:]
            try:
                pos = int(row['DoE position'])
                part = row['Part number']
            except:
                continue
            if pos not in levels:
                levels[pos] = []
            levels[pos].append( part )
        df = pd.read_csv(arg.plasmidFile)
        des = []
        for i in np.arange(df.shape[0]):
            if i > int(arg.n):
                break
            const = []
            newcol = []
            for col in df.columns:
                try:
                    posi = int( col.split('.')[0] )
                    rol = col.split('.')[1]
                except:
                    break
                if posi in levels and len(levels[posi]) > 1:
                    const.append( "L{}".format( levels[posi].index(df.loc[i,col]) +1 ) )
                    newcol.append( rol+str(posi) )
            des.append( const )
        dfdes = pd.DataFrame(des, columns=newcol)
        dfdes.to_csv(os.path.join(os.path.dirname(arg.plasmidFile),
                        re.sub('.csv','',os.path.basename(arg.plasmidFile))+'_newDoE.csv'),
                        index=False)
            
            
        
    
