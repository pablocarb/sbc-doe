'''
mapSamples (c) University of Manchester 2018

mapSamples is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Map samples into their DoEs ids (and generate factors table)
@usage: 
    - Analyse designs: mapSamples.py dtsFolder -outFolder OUTFOLDER -iceEntries ICEENTRIES -designsFolder DESIGNSFOLDER
    - Generate only the table: mapSamples.py dtsFolder -utFolder outFolder -iceEntries iceEntriesFile -onlyTable
'''
import os, re, argparse, csv, glob
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from itertools import product, permutations
from synbiochem.utils import ice_utils
import sys

LIVEBUTTON = False

def fullImport(arg):
    """ Import only for design evaluation """
    global evalDesignInfo
    if not arg.onlyTable:
        sys.path.append(os.path.join(os.getenv('CODE'),'sbml2doe'))
        from stress_test_scripts.evaluateSBCDesign import evalDesignInfo
        

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
    parser.add_argument('-onlyTable', action='store_true',
                        help='Generate only table')
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

def patch31():
    plmap = {
             'SBC007160': 'SBCDE00031_PL01', 
             'SBC007162': 'SBCDE00031_PL02',
             'SBC007164': 'SBCDE00031_PL03',
             'SBC007166': 'SBCDE00031_PL04',
             'SBC007168': 'SBCDE00031_PL05',
             'SBC007170': 'SBCDE00031_PL06',
             'SBC007172': 'SBCDE00031_PL07',
             'SBC007174': 'SBCDE00031_PL08',
             'SBC007176': 'SBCDE00031_PL09',
             'SBC007178': 'SBCDE00031_PL10',
             'SBC007180': 'SBCDE00031_PL11',
             'SBC007182': 'SBCDE00031_PL12',
             'SBC010238': 'SBCDE00049_PL01'
             }
    return plmap


def mapPlasmids(df, arg, columns=['Strain ID','Plasmid ID'], mapColumn='Plasmid Name', iceurl="https://ice.synbiochem.co.uk"):
    """ Map plasmids into their names in ICE, they should contain the DoE plasmid name """
    plmap = {}
    if arg.iceEntries is not None and os.path.exists( arg.iceEntries ):
        icedf = pd.read_csv( arg.iceEntries )
        # Patch depending if the entries file is retrieved from json or from the website
        if 'partId' in icedf.columns:
            partIdix = 'partId'
            nameix = 'name'
        else:
            partIdix = 'Part ID'
            nameix = 'Name'
        for i in icedf.index:
            plmap[ icedf.loc[i,partIdix] ] = icedf.loc[i,nameix]
    else:
        ids = set()
        for column in columns:
            try:
                for x in df.loc[:,column]:
                    sbcid = filterSBCid(x)
                    if sbcid is not None:
                        ids.add( sbcid )
            except:
                continue
        client = ice_utils.ICEClient(iceurl,os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
        plmap = {}
        for sbcid in ids:
            try:
                entry = client.get_ice_entry( sbcid )
                name = entry.get_name()
                plmap[ sbcid ] = name
            except:
                continue
    patch = patch31()
    for x in patch:
        plmap[x] = patch[x]
    for column in columns:
        for i in range(0, df.shape[0]):
            try:
                x = df.loc[i,column]
                sbcid = filterSBCid(x)
                if sbcid is not None and sbcid in plmap:
                    df.loc[i, mapColumn] = plmap[ sbcid ]
            except:
                continue
    return df

def outputSamples(df, outputFolder):
    outfile = os.path.join(outputFolder, 'samples.csv')
    df.to_csv( outfile )
    

def getDesignSpec(doeinfo):
   # Look for allowed positions for each part
    constraints = {}
    movpos = set()
    movpart = set()
    for i in doeinfo.index:
        part = doeinfo.loc[i,'Part number']
        pos = doeinfo.loc[i,'DoE position']
        ptype = doeinfo.loc[i,'DoE designation']
        perm = doeinfo.loc[i,'DoE permutation class']
        if part == '-':
            part = 'None'
        try:
            pos = int(pos)
        except:
            continue
        if perm == 1.0:
            movpos.add( pos )
            movpart.add( part )
        if pos not in constraints:
            constraints[pos] = set()
        constraints[pos].add( part )
    for pos in movpos:
        for part in movpart:
            constraints[pos].add( part )
    # const2 are the factor constraints, i.e which factors are possible at each position
    # map2 are possible values that can take the factor with no constraints 
    # (to be used only if the predicted DoE has not introduced a constraint )
    const2 = {}
    map2 = {}
    cols = {}
    combi = {}
    construct = {}
    for i in doeinfo.index:
        part = doeinfo.loc[i,'Part number']
        pos = doeinfo.loc[i,'DoE position']
        ptype = doeinfo.loc[i,'DoE designation']
        perm = doeinfo.loc[i,'DoE permutation class']
        if part == '-':
            part = 'None'
        try:
            pos = int(pos)
        except:
            continue
        cols[pos] = ptype
        if pos not in combi:
            combi[pos] = set()
        combi[pos].add(part)
        partclass = set()
        if ptype == 'origin':
            pid = 'pl'
            partclass.add( pid )
            construct[pos] = pid
        elif ptype == 'resistance':
            pid = 'res'
            partclass.add( pid )
            construct[pos] = pid
        elif ptype == 'promoter':
            pnext = pos + 1
            # if the gene next is mobile
            # then add
            if pnext in movpos:
                for p in movpos:
                    pid = 'p_'+'g'+str(p)
                    partclass.add( pid )
            else:
                pid = 'p_'+'g'+str(pnext)
                partclass.add( pid  )
            pid = 'p_'+'g'+str(pnext)   
            construct[pos] = pid
    
        # Not clear how to deal with multple variants: to do
        elif ptype == 'gene':
            if pos in movpos:
                for p in movpos:
                    pid =  'g_'+'g'+str(p)
                    partclass.add( pid )
            else:
                pid = 'g_'+'g'+str(pos)
                partclass.add( pid )
            pid = 'g_'+'g'+str(pos)
            construct[pos] = pid
        const2[pos] = partclass
        for p in partclass:
            if p not in map2:
                map2[p] = set()
            map2[p].add(part)
    construct = [ construct[x] for x in sorted(construct)]
    return constraints, const2, map2, cols, movpos, combi, construct
    
def bestPlasmids(ndata, doeinfo):
    """ Start with best predicted combinations
    convert into possible plasmids """
    # Look for allowed positions for each part
    constraints, const2, map2, cols, movpos, combi, construct = getDesignSpec(doeinfo)
    nplas = []
    for index, row in ndata.iterrows():
        w = []
        watch = set()
        pwatch = []
        for pos in sorted( const2 ):
            v = []
            # Loop through all known instances at the position
            # Sanity check, i.e. whether:
            # - They are allowed according to the constraints
            # - They have not been used yet
            for x in sorted(const2[pos]):
                if x in row:
                    if row[x] in constraints[pos] and row[x] not in watch:
                        if x.startswith('g_'):
                            promo = re.sub('g_','p_',x)
                            if promo in row and promo != pwatch[-1]:
                                continue
                        v.append( row[x] )  
                        watch.add(row[x])
                        pwatch.append(x)
                else:
                    for y in map2[x]:
                        if y not in watch:
                            v.append(y)
                            watch.add(y)
                            pwatch.append(x)
            w.append(v)
        for h in product( *w ):
            nplas.append(h+(row['pred'],row['Exp'],row['Design']))
    coln = []
    px = 1
    for z in sorted(cols):
        coln.append( "{}.{}".format(px, cols[z]))
        px += 1
    coln += [ 'pred', 'Exp', 'Design' ]
    return pd.DataFrame( nplas, columns=coln )

def spatArrang(construct, movpos):
    consts = []
    mobiles = sorted(movpos)
    perm = [x for x in permutations( mobiles )]
    for v in perm:
        c1 = np.arange(len(construct))
        for y in np.arange(len(v)):
            c1[mobiles[y]-1] = v[y]-1
        consts.append(c1)
    return consts

def getConstruct(construct, z):
    """ After rearrangement of genes, promoters are relabelled if they are present """
    geneArr = np.array(construct)[z]
    proms = set()
    for x in geneArr:
        if x.startswith('p_'):
            proms.add(x)
    for i in np.arange(len(geneArr)):
        if geneArr[i].startswith('g_'):
            pg = re.sub('g_', 'p_', geneArr[i])
            # Make sure that the promoter exist in the construct
            if pg in proms and geneArr[i-1].startswith('p_'):
                geneArr[i-1] = re.sub('g_', 'p_', geneArr[i])
    return geneArr


def filterUnseen(df, ndata, descol):
    """ Remove cases containing categorical variables not used in the training set.
    This happens often with promoters in front of mobile genes: not all of them occur """
    pathol = []
    for x in descol:
        s1 = set(df[x].unique())
        s2 = set(ndata[x].unique())
        for unk in s2-s1:
            pathol.append( ndata[x] != unk )
    if len(pathol) > 0:
        ix = np.logical_and.reduce( pathol )
        ndata = ndata.loc[ix,:]
    return ndata

def mapCombi(nc, arr, descol):
    """ Map combinations into the future vector """
    feat = [0 for i in np.arange(len(descol))]
    for i in np.arange(len(arr)):
        v = arr[i]
        val = nc[i]
        if v in descol:
            feat[ descol.index(v) ] = val
            if val == 'None':
                g1 = re.sub('g_','', arr[i-1])
                g2 = re.sub('g_','', arr[i+1])
                v1 = '_'.join([g1,g2])
                if v1 in descol:
                    feat[ descol.index(v1) ] = 1
            elif arr[i].startswith('g_'):
                try:
                    assert arr[i+1].startswith('g_')
                except:
                    continue
                g1 = re.sub('g_','', arr[i])
                g2 = re.sub('g_','', arr[i+1])
                v1 = '_'.join([g1,g2])
                if v1 in descol:
                    feat[ descol.index(v1) ] = 1
                
    return feat
    

def bestCombinations(df, res, doeinfo=None, version=2):
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
    VERSION = 2
    # If multiple positional combinations are possible
    if len(pos) > 0:
        if version == 1:
            # 1st version tries to predict known combinations
            # However, it does not work very well
            ucol = df.loc[:,pos].drop_duplicates()
            # For each spatial arrangement
            for x in ucol.index:
                # Select all rows containing the arrangement
                ix = None
                for y in ucol.columns:
                    cond = df[y] == ucol.loc[x,y]
                    if ix is None:
                        ix = cond
                    else:
                        ix = np.logical_and(ix, cond)
                ugroup.append(ix)
                # Store all instances for the spatial arrangement
                combo = []
                for z in descol:
                    combo.append( df.loc[ix,z].unique() )
                ucomb.append( combo )
                combi.extend( [w for w in product( *combo )] )
        elif VERSION == 2:
            # To do: this is called twice
            constraints, const2, map2, cols, movpos, comb, construct = getDesignSpec(doeinfo)
            # Generate the full combinations without spatial rearrangements
            combo = []
            for pos in sorted(comb):
                combo.append( comb[pos] )
            combi = [w for w in product( *combo )]
            # To DO: ADD SPATIAL ARRANGEMENTS
            consts =  spatArrang(construct, movpos)
            ncombi = []
            for w in combi:
                for z in consts:
                    nc = np.array(w)[z]
                    arr = getConstruct(construct, z)
                    feat = mapCombi(nc,arr,descol)
                    ncombi.append( feat )
            combi = ncombi
    else:
            combo = []
            for z in descol:
                combo.append( df.loc[:,z].unique() )
            ucomb.append( combo )
            combi.extend( [w for w in product( *combo )] )
    ndata = pd.DataFrame(combi, columns=descol)
    ndata = filterUnseen(df, ndata, descol)
    ndata['pred'] = res.predict( ndata )
    ndata = ndata.sort_values(by='pred', ascending=False)
    return ndata

def title(text, level=2):
    return '<h{}>{}</h{}>'.format(str(level), text, str(level))

def livetitle(text, nsection, level=2):
    string = '<h{}><a  href="#" onclick="'+"toggle_visibility('sec{}');"+'" title="Show/Hide">+ </a>{}</h{}>'
    return string.format(str(level), nsection, text, str(level))

def startSection( nsection, sectitle ):
    if LIVEBUTTON:
        ht = livetitle( sectitle, nsection )
    else:
        ht = title( str(nsection+1)+'. '+sectitle )        
    ht += '<div id="sec{}">'.format( nsection )
    return ht

def statsHTML(outfile, info, desid, target, doeinfo, ndata, nplasm, desinfo=[]):
    """ Output an HTML file with the analysis """
    """ target: targets[t] """
    htm = '<link rel="stylesheet" href="style.css">'
    htm += '<script src="script.js"></script>'
    htm += title('Predictive Analytics', 1)
    htm += '<div><b>Design: </b>'+desid+'</div>'
    htm += '<div><b>Target: </b>'+target+'</div>'
    sections = ['Model fitting:', 'Contrast and regression effects:', 'Model diagnostics:']
    nTable = 0
    for line in info.as_html().split('\n'):
        if line.startswith('<table'):
#            htm += title(sections[nTable])
            htm += startSection( nTable, sections[nTable] ) 
            nTable += 1
        htm += line
        if line.startswith('</table'):
            htm += '</div>'
    if doeinfo is not None:
        ix = []
        for i in doeinfo.index:
            try:
                int( doeinfo.iloc[i,0] )
                ix.append( i )
            except:
                continue
#        htm += title('DoE specifications:')
        htm += startSection( nTable, 'DoE specifications:' )
        nTable += 1
        htm += doeinfo.iloc[ix,0:7].to_html(index=False)
        htm += '</div>'
#    htm += title('Predicted best combinations:')
    htm += startSection( nTable, 'Predicted best combinations:' )
    nTable += 1
    htm += ndata.to_html(index=False)
    htm += '</div>'
    if nplasm is not None:
#        htm += title('Predicted best constructs:')
        htm += startSection( nTable, 'Predicted best constructs:' )
        nTable += 1
        htm += nplasm.to_html(index=False)
        htm += '</div>'
    if len(desinfo) > 0:
        tab1 = desinfo[0]
        htm += startSection( nTable, 'Design evaluation:' )
        nTable += 1
        htm += tab1.to_html(index=False)
        htm += '</div>'
    if len(desinfo) > 2 and desinfo[2] is not None:
        tab1 = desinfo[2]
        htm += startSection( nTable, 'Factor power analysis:' )
        nTable += 1
        htm += tab1.to_html(index=False)
        htm += '</div>'
    if len(desinfo) > 1 and desinfo[1] is not None:
        tab1 = desinfo[1]
        htm += startSection( nTable, 'Samples relative prediction variance:' )
        nTable += 1
        htm += tab1.to_html(index=False)
        htm += '</div>'
        
    with open(outfile, 'w') as h:
        h.write(htm)

def stats(df, desid, doeinfo=None, outputFolder='/mnt/SBC1/data/Biomaterials/learn', desinfo=[]):
    """ Perform some statistical analysis of the factors 
        - df: DataFrame containing the factors and the response;
        - doeinfo: DataFrame generated from the info file;
        - outputFolder
    """
    targets = {}
    targetsList = []
    for j in np.arange(len(df.columns)):
        x = df.columns[j]
        if x.startswith('Target') and x.endswith('Conc'):
            if len( df[x].unique() ) > 1:
                targets[ x ] = None
                for i in df.index:
                    try:
                        if np.isnan( df.loc[i,df.columns[j-1]] ):
                            continue
                    except:
                        pass
                    targets[ x ] = df.loc[i,df.columns[j-1]]
                    targetsList.append( x )
                    break
            
    factors = []
    for x in np.flip( df.columns ):
        if x == 'Design':
            break
        factors.append( x )
    factors.reverse()
    ftargets = {}
    for i in np.arange(len(targetsList)):
        t = targetsList[i]
        # Skip columns with no named target
        if targets[t] is None:
            continue
        formula = "Q('{}') ~".format( t )
        terms = []
        for f in factors:
            terms.append( "Q('{}')".format( f ) )
        formula += ' + '.join(terms)
        ols = smf.ols( formula=formula, data=df)
        res = ols.fit()
        info = res.summary()
        ndata = bestCombinations(df, res, doeinfo)
        match = df.drop_duplicates('Design')
        # Merge with predictions, preserve left index
        vv = pd.merge( ndata, match, how='inner', on=list(ndata.columns)[0:-1], validate='one_to_one', right_index=True) 
        # Calculate experimental means (o median, etc)
        mm = df.groupby('Design').mean()
        for k in vv.index:
            ndata.loc[k,'Exp'] = float( mm.loc[ mm.index == vv.loc[k,'Design'], t] )
            ndata.loc[k,'Design'] = vv.loc[k,'Design']
        if doeinfo is not None:
            nplasm = bestPlasmids(ndata, doeinfo)
        else:
            nplasm = None
        # Keep 1 row per design

        outfile = os.path.join( outputFolder, desid+'_summary_'+str(i)+'.csv' )
        cv = 'Design: , {}\n'.format(desid)
        cv +='Target: , {}\n'.format(targets[t])
        cv += info.as_csv()
        with open(outfile, 'w') as h:
            h.write( cv )
        outinfo = desid+'_summary_'+str(i)
        outname =  outinfo+'.html' 
        outfile = os.path.join( outputFolder,outname )
        ftargets[t] = outname
        statsHTML(outfile, info, desid, targets[t], doeinfo, ndata, nplasm, desinfo=desinfo)
        ndataf = os.path.join( outputFolder,outinfo+'_constructs.csv')
        ndata['Target'] = targets[t]
        ndata.to_csv(ndataf,index=False)
        nplasm['Target'] = targets[t]
        nplasmf = os.path.join( outputFolder,outinfo+'_plasmids.csv')
        nplasm.to_csv(nplasmf,index=False)
    return targetsList, targets, ftargets

def outputFactors(df, designsFolder='/mnt/syno/shared/Designs',
                  outputFolder='/mnt/SBC1/data/Biomaterials/learn', mapColumn='Plasmid Name',
                  plateColumn='Plate ID', info=None ):
    """ Output the file with the factors and call the learn routine """
    """ Info: add extra info like file name in order to identify better the file """
    outcome = []
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
        # Read design 
        jfile = os.path.join( designsFolder, des,'Design', 'DoE_'+re.sub('SBCDE', 'SBC', des)+'.dat')
        if os.path.exists(jfile):
            tab1, tab2, tab3 = evalDesignInfo(jfile)
            desinfo = [tab1, tab2, tab3]
        else:
            desinfo = []
        # Add "independent" factors
        facdict = {}
        for i in fcdf.index:
            facdict[fcdf.loc[i,'Design']] = i
        for i in desId[des]:
            design = df.loc[i,mapColumn]
            if design in facdict:
                rows.append( np.hstack( [df.loc[i,:],fcdf.loc[facdict[design],:]] ) )
        # Perform stats, output results
        desidrnd = str( np.random.randint(1000) )
        if len(rows) > 0:
            fulldf = pd.DataFrame( rows, columns=np.hstack( [df.columns, fcdf.columns] ) )
            if info is None:
                desi = des+'_'+'UNKNOWN'+str( desidrnd )
                if des in plateId:                   
                    try:
                        desi = des+'_'+plateId[des]
                    except: 
                        pass
            else:
                desi = des+'_'+info
            
            targetsList, targets, ftargets = stats( fulldf, desi, doeinfo, outputFolder, desinfo=desinfo )
            if len(targetsList) > 0:
                print( targetsList )
                outcsv = os.path.join(outputFolder, desi+'_learn.csv')
                fulldf.to_csv( outcsv )
                for t in targetsList:
                    outcome.append( (des,targets[t],ftargets[t]) )
    return outcome
    
def readDesign( dfile, des={} ):
    with open(dfile) as h:
        for line in h:
            row = line.rstrip().split('\t')
            des[ row[0] ] = row[1:]

def addDesignColumns( df, designsFolder='/mnt/syno/shared/Designs', ext='.j0', mapColumns=['Plasmid Name','Plasmid ID']):
    """ Create additional columns with the design combinations """
    """ TO DO: transform positional parameters into: promoters in front of genes, gene variants, and pairings  """
    for mapColumn in mapColumns:
        designs = set()
        for plasmid in df[mapColumn].unique():
            try:
                if plasmid.startswith('SBCDE'):
                    desname = plasmid.split('_')[0]
                    desno = int( re.sub('SBCDE','',desname) )
                    desname = "SBCDE%05d" % (desno,)
                    designs.add( desname )
            except:
                continue
        des = {}
        for desi in designs:
            dfile = os.path.join( designsFolder, desi, 'Design', desi+ext )
            readDesign( dfile, des )
        for i in df.index:
            plasmid = df.loc[i,mapColumn]
            if plasmid in des:
                for j in range(0, len( des[plasmid] ) ):
                    df.loc[i,'pos'+str(j)] = des[plasmid][j]
 
def makeSummary(outcome):
    htm = '<link rel="stylesheet" href="style.css">'
    summary = pd.DataFrame(columns=['Design','Target','DTS','Link'])    
    for target in outcome:
        row = dict( zip(summary.columns, [ target[1],target[2],target[0],'<a href="{}" target="_blank">{}</a>'.format(target[3],target[3]) ]))
        summary = summary.append( row, ignore_index=True )
    summary.sort_values( by=['Design','Target','DTS'], inplace=True )
    pd.set_option('display.max_colwidth',-1)
    htm += summary.to_html(escape=False,index=False)
    summary.to_excel(os.path.join(arg.outFolder,'index.xlsx'))
    with open(os.path.join(arg.outFolder,'index.html'), 'w') as h:
        h.write(htm)

def readPlate(dts,sheet='Media'):
    val = []
    d1 = pd.read_excel(dts, sheet)
    for plate in np.arange(0,5):
        for i in np.arange(0,8):
            row = i + plate*11
            for j in np.arange(0,12):
                col = j+1
                try:
                    val.append(d1.iloc[row,col])
                except:
                    continue
    return val

def addSupplInfo(df, dts, noEmptyTargets=False):
    """ Add additional columns and filter non-relevant columns """
    if noEmptyTargets:
        drop = set()
        for col in df.columns:
            if col.startswith('Target'):
                tarid = ' '.join(col.split(' ')[0:2])
                if col.endswith('Name'):
                    try:
                        if np.isnan( df.loc[0,col] ):
                            drop.add(tarid)
                    except:
                        continue
        dropList = []
        for col in df.columns:
            if col.startswith('Target'):
                tarid = ' '.join(col.split(' ')[0:2])
                if tarid in drop:
                    dropList.append(col)
        df = df.drop(dropList, axis=1)
    for col in ['Strain','Plasmid','Media','Treatment']:
        if col not in df.columns:
            val = readPlate(dts,sheet='Media')
            df[col] = val
    for row in df.index:
        pid = df.loc[row,'Plasmid Name']
        try:
            if pid.startswith('SBCDE'):
                sbcd, pln = pid.split('_')
                df.loc[row,'Design ID'] = sbcd
        except:
            df.loc[row,'Design ID'] = None
    for row in df.index:
        pid = df.loc[row,'Plasmid ID']
        try:
            if np.isnan(pid):
                df.loc[row,'Plasmid ID'] = df.loc[row,'Strain ID']
        except:
            continue
    df['dts'] = os.path.basename(dts)  
    cols = ['dts','Plate Number', 'Well', 'Plate ID', 'Description', 'Strain ID',
            'Strain Name', 'Host ID', 'Host Name', 'Plasmid ID', 'Plasmid Name',
            'Design ID', 'Media', 'Treatment', 'OD Induction', 'OD Harvest', 
            'Date Created', 'Creator']
    for col in df.columns:
        if col.startswith('Target'):
            cols.append(col)
    for missing in set(cols) - set(df.columns):
        df[missing] = None
    df = df[cols]
    nodata = []
    for row in df.index:
        empty = True
        for col in df.columns:
            if col.endswith('Conc'):
                try:
                    if not np.isnan( df.loc[row,col] ):
                        empty = False
                except:
                    empty = False
        if empty:
            nodata.append( row )
    df = df.drop( nodata, axis=0 )
    return df
    

if __name__ == '__main__':
    parser = arguments()    
    arg = parser.parse_args()
    fullImport(arg)
    big = None
    outcome = []
    for dirName in glob.glob( os.path.join(arg.dts,'*') ):
        for dts in glob.glob( os.path.join(dirName, '*') ):
            if not (dts.endswith('.xlsm') or dts.endswith('.XLSM')):
                continue
            if dts.lower().endswith( 'xlsm') and not os.path.basename(dts).startswith('~'):# and len(dts.split('_')) == 1:
                print( dts )
                if not dts.endswith('190814QQQMAN07.XLSM'):
                    continue
                # read a DataFrame with the samples
                df = samples( os.path.join(dirName, dts) )
                noEmpty = False
                for col in df.columns:
                    if col.endswith('Conc'):
                        noEmpty |= not np.all(np.isnan(df['Target 1 Conc']))
                if not noEmpty:
                    continue
                print( 'Reading data' )
                # Map plasmids into their names in ICE, they should contain the DoE plasmid name
                df = mapPlasmids( df, arg )
                df = addSupplInfo(df, os.path.join(dirName, dts))
                if big is None:
                    big = df
                else: 
                    big = big.append(df, ignore_index=True)
                if arg.onlyTable:
                    continue
                # Create additional columns with the design combinations
                addDesignColumns( df, designsFolder=arg.designsFolder )
            #    outputSamples( df, outputFolder=arg.outFolder )
                # Output some statistics about the factors
                vals = outputFactors( df, designsFolder=arg.designsFolder, outputFolder=arg.outFolder ) 
                location = dirName.split(arg.dts)[1]
                for x in vals:
                    outcome.append( (os.path.join(location,os.path.basename(dts)),)+x )
    if not arg.onlyTable:
        makeSummary(outcome)
    if big is not None:
        big.to_excel(os.path.join(arg.outFolder,'bigtable.xlsx'), index=False)
