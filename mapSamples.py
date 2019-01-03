'''
mapSamples (c) University of Manchester 2018

mapSamples is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Map samples into their DoEs ids (and generate factors table)
@usage: mapSamples.py dts -outFolder OUTFOLDER -iceEntries ICEENTRIES -designsFolder DESIGNSFOLDER
'''
import os, re, argparse, csv
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from itertools import product
from synbiochem.utils import ice_utils

def arguments():
    parser = argparse.ArgumentParser(description='MapSamples: Learn design rules. Pablo Carbonell, SYNBIOCHEM, 2018')
    parser.add_argument('dts', 
                        help='Data tracking sheet folder')
    parser.add_argument('-outFolder',  default='/mnt/SBC1/data/Biomaterials/learn',
                        help='Output folder ')
    parser.add_argument('-iceEntries',  default=None,
                        help='ICE entries csv file (instead of accessing client)')
    parser.add_argument('-designsFolder', default='/mnt/syno/shared/Designs',
                        help='Folder that contains the design templates')
    return parser

def samples(dts, sheetName='Samples'):
    try:
        df = pd.read_excel(dts, sheetName)
    except:
        df = None
    return df

def filterSBCid(x):
    """ Try to figure out the SBC id """
    try:
        iceid = int(x)
        sbcid = 'SBC'+"%06d" % (iceid,)
    except:
        try:
            if x.startswith('SBC') and len(x) >= 9:
                sbcid = x[0:9]
            else:
                sbcid = None
        except:
            sbcid = None
    return sbcid

def mapPlasmids(df, arg, column='Strain ID', mapColumn='Plasmid Name', iceurl="https://ice.synbiochem.co.uk"):
    """ Map plasmids into their names in ICE, they should contain the DoE plasmid name """
    plmap = {}
    if arg.iceEntries is not None and os.path.exists( arg.iceEntries ):
        icedf = pd.read_csv( arg.iceEntries )
        for i in icedf.index:
            plmap[ icedf.loc[i,'Part ID'] ] = icedf.loc[i,'Name']
    else:
        ids = set()
        for x in df.loc[:,column]:
            sbcid = filterSBCid(x)
            if sbcid is not None:
                ids.add( sbcid )
        client = ice_utils.ICEClient(iceurl,os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
        plmap = {}
        for sbcid in ids:
            try:
                entry = client.get_ice_entry( sbcid )
                name = entry.get_name()
                plmap[ sbcid ] = name
            except:
                continue
    for i in range(0, df.shape[0]):
        x = df.loc[i,column]
        sbcid = filterSBCid(x)
        if sbcid is not None and sbcid in plmap:
            df.loc[i, mapColumn] = plmap[ sbcid ]
        else:
            df.loc[i, mapColumn] = 'None'
    return df

def outputSamples(df, outputFolder):
    outfile = os.path.join(outputFolder, 'samples.csv')
    df.to_csv( outfile )

def bestCombinations(df, res):
    """ Predict best allowed combinations (experimental)
    Based solely on the values present on the library
    (no suggestion for introducing new changes)
    """
    descol = [i for i in df.columns]
    descol.reverse()
    
    for i in np.arange(len(descol)):
        if descol[i] == 'Design':
            break
    descol = descol[:i]
    descol.reverse()
    pos = []
    for d in descol:
        try:
            a,b = d.split('_')
        except:
            continue
        if re.match('g\d+',a) and re.match('g\d+',b):
            pos.append(d)
    # unique combinations
    combi = []
    ugroup = []
    ucomb = []
    # If multiple positional combinations are possible
    if len(pos) > 0:
        ucol = df.loc[:,pos].drop_duplicates()
        for x in ucol.index:
            ix = None
            for y in ucol.columns:
                cond = df[y] == ucol.loc[x,y]
                if ix is None:
                    ix = cond
                else:
                    ix = np.logical_and(ix, cond)
            ugroup.append(ix)
            combo = []
            for z in descol:
                combo.append( df.loc[ix,z].unique() )
            ucomb.append( combo )
            combi.extend( [w for w in product( *combo )] )
    else:
            combo = []
            for z in descol:
                combo.append( df.loc[:,z].unique() )
            ucomb.append( combo )
            combi.extend( [w for w in product( *combo )] )
    ndata = pd.DataFrame(combi, columns=descol)
    ndata['pred'] = res.predict( ndata )
    ndata = ndata.sort_values(by='pred', ascending=False)
    return ndata

def title(text, level=2):
    return '<h{}>{}</h{}>'.format(str(level), text, str(level))

def stats(df, desid, doeinfo=None, outputFolder='/mnt/SBC1/data/Biomaterials/learn'):
    """ Perform some statistical analysis of the factors 
        - df: DataFrame containing the factors and the response;
        - doeinfo: DataFrame generated from the info file;
        - outputFolder
    """
    targets = {}
    for j in np.arange(len(df.columns)):
        x = df.columns[j]
        if x.startswith('Target') and x.endswith('Conc'):
            if len( df[x].unique() ) > 1:
                targets[ x ] = 'Unknown compound'
                for i in df.index:
                    try:
                        if np.isnan( df.loc[i,df.columns[j-1]] ):
                            continue
                    except:
                        pass
                    targets[ x ] = df.loc[i,df.columns[j-1]]
                    break
            
    factors = []
    for x in np.flip( df.columns ):
        if x == 'Design':
            break
        factors.append( x )
    factors.reverse()
    targetsList = sorted( targets )
    for i in np.arange(len(targetsList)):
        t = targetsList[i]
        formula = "Q('{}') ~".format( t )
        terms = []
        for f in factors:
            terms.append( "Q('{}')".format( f ) )
        formula += ' + '.join(terms)
        ols = smf.ols( formula=formula, data=df)
        res = ols.fit()
        info = res.summary()
        ndata = bestCombinations(df, res)
        outfile = os.path.join( outputFolder, desid+'_summary_'+str(i)+'.csv' )
        cv = 'Design: , {}\n'.format(desid)
        cv +='Target: , {}\n'.format(targets[t])
        cv += info.as_csv()
        with open(outfile, 'w') as h:
            h.write( cv )
        outfile = os.path.join( outputFolder, desid+'_summary_'+str(i)+'.html' )
        htm = '<link rel="stylesheet" href="style.css">'
        htm += title('Predictive Analytics', 1)
        htm += '<div><b>Design: </b>'+desid+'</div>'
        htm += '<div><b>Target: </b>'+targets[t]+'</div>'
        sections = ['Model fitting:', 'Contrast and regression effects:', 'Model diagnostics:']
        nTable = 0
        for line in info.as_html().split('\n'):
            if line.startswith('<table'):
                htm += title(sections[nTable])
                nTable += 1
            htm += line
        if doeinfo is not None:
            ix = []
            for i in doeinfo.index:
                try:
                    int( doeinfo.iloc[i,0] )
                    ix.append( i )
                except:
                    continue
            htm += title('DoE specifications:')
            htm += doeinfo.iloc[ix,0:7].to_html(index=False)
        htm += title('Predicted best combinations:')
        htm += ndata.to_html(index=False)
        with open(outfile, 'w') as h:
            h.write(htm)

def outputFactors(df, designsFolder='/mnt/syno/shared/Designs',
                  outputFolder='/mnt/SBC1/data/Biomaterials/learn', mapColumn='Plasmid Name',
                  plateColumn='Plate ID'):
    plateId = {}
    desId = {}
    # Loop through the dts, retrieve design + plasmid id from the mapColumn, 
    # retrieve plateid from plateColumn
    for i in df.index:
        plid = df.loc[i,mapColumn]
        if type(plid) == str and plid.startswith('SBCDE'):
            try:
                sbcid, plid = plid.split('_')
            except:
                continue
            if sbcid not in desId:
                desId[sbcid] = []
            desId[sbcid].append( i )
            plateId[sbcid] = df.loc[i,plateColumn]
    # For each design, read factors, read doe, add "independent" factors, perform stats
    for des in desId:
        # Read factors
        rows = []
        factorFile = os.path.join(designsFolder, des, 'Design', des+'_factors.csv')
        if os.path.exists(factorFile):
            fcdf = pd.read_csv(factorFile)
        else:
            continue
        # Read doe
        doeFile = os.path.join(designsFolder, des, 'Design', 'DoE_'+re.sub('SBCDE', 'SBC', des)+'.xlsx')
        if os.path.exists(doeFile):
            doeinfo = pd.read_excel( doeFile )
        else:
            doeinfo = None
        # Add "independent" factors
        facdict = {}
        for i in fcdf.index:
            facdict[fcdf.loc[i,'Design']] = i
        for i in desId[des]:
            design = df.loc[i,mapColumn]
            if design in facdict:
                rows.append( np.hstack( [df.loc[i,:],fcdf.loc[facdict[design],:]] ) )
        # Perform stats, output results
        if len(rows) > 0:
            fulldf = pd.DataFrame( rows, columns=np.hstack( [df.columns, fcdf.columns] ) )
            desi = des+'_'+plateId[des]
            outcsv = os.path.join(outputFolder, desi+'_learn.csv')
            fulldf.to_csv( outcsv )
            stats( fulldf, desi, doeinfo, outputFolder )

def readDesign( dfile, des={} ):
    with open(dfile) as h:
        for line in h:
            row = line.rstrip().split('\t')
            des[ row[0] ] = row[1:]

def addDesignColumns( df, designsFolder='/mnt/syno/shared/Designs', ext='.j0', mapColumn='Plasmid Name'):
    """ Create additional columns with the design combinations """
    """ TO DO: transform positional parameters into: promoters in front of genes, gene variants, and pairings  """
    designs = set()
    import pdb
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
    for dirName, subdirList, fileList in os.walk( arg.dts ):
        for dts in fileList:
            if dts.lower().endswith( 'xlsm') and not dts.startswith('~'):# and len(dts.split('_')) == 1:
                print( dts )
                # read a DataFrame with the samples
                df = samples( os.path.join(dirName, dts) )
                # Map plasmids into their names in ICE, they should contain the DoE plasmid name 
                df = mapPlasmids( df, arg )
                # Create additional columns with the design combinations
                addDesignColumns( df, designsFolder=arg.designsFolder )
            #    outputSamples( df, outputFolder=arg.outFolder )
                # Output some statistics about the factors
                outputFactors( df, designsFolder=arg.designsFolder, outputFolder=arg.outFolder )
