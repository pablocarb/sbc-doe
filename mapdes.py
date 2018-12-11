'''
mapdes (c) University of Manchester 2018

mapdes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell
@description: Map design of experiments to SYNBIOCHEM registered designs by querying ICE
'''
import os, re, argparse, csv
from synbiochem.utils import ice_utils


def arguments():
    parser = argparse.ArgumentParser(description='MapDes. Pablo Carbonell, SYNBIOCHEM, 2018')
    parser.add_argument('designFile', 
                        help='DOE design file')
    parser.add_argument('outFile', 
                        help='DOE mapping output file')
    return parser

def loop():
    outbase = '/mnt/SBC1/data/Biomaterials/mappings'
    base = '/mnt/syno/shared/Designs'
    client = ice_utils.ICEClient("https://ice.synbiochem.co.uk",os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
    for des in range(31, 50):
        name = 'SBCDE000'+str(des)
        desfile = os.path.join( base, name, 'Design', name+'.txt' )
        if not os.path.exists( desfile ):
            continue
        outfile = os.path.join( outbase, name+'_map.csv' )
        grabInfo( client, desfile, outfile )

def grabInfo( client, infile, outfile ):
    folder = os.path.dirname( infile )
    with open( infile ) as h, open( outfile, 'w') as w:
        cw = csv.writer(w)
        cw.writerow( ( 'plasmid', 'partId', 'designId', 'fullName', 'shortDescription' ) )
        for row in h:
            m = row.rstrip().split('\t')
            plasmid = m[0]
            result = client.search( plasmid )
            for r in result['results']:
                if r['entryInfo']['type'] == 'PLASMID':
                    partId = r['entryInfo']['partId']
                    short = r['entryInfo']['shortDescription']
                    entry = client.get_ice_entry( partId )
                    param = entry.get_parameters()
                    try:
                        designId = param['Design id']
                    except:
                        designId = None
                    try:
                        fullName = param['Full name']
                    except:
                        fullName = None
                    row =  (plasmid, partId, designId, fullName, short)
                    print( row )
                    cw.writerow( row )
                

if __name__ == '__main__':
    parser = arguments()
    arg = parser.parse_args()
    client = ice_utils.ICEClient("https://ice.synbiochem.co.uk",os.environ['ICE_USERNAME'], os.environ['ICE_PASSWORD'])
#    mapfile = os.path.join( folder, re.sub('\.txt', '.map', os.path.basename( arg.designFile ) ) )
    grabInfo( client, arg.designFile, arg.outFile )
