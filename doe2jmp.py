'''
Doe2JMP (c) University of Manchester 2018

Doe2JMP is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Prepare files for JMP from the DoE sheet
'''
from doe import read_excel
import argparse
import os, re, sys
import numpy as np

def makeDoeScript(fact, outfile, size, seed=None, starts=1040, executable=False, makeTable=False):
    """ Full Doe Script in JMP """
    doe = []
    if executable:
        doe.append( '//!' )
    doe.append( 'DOE(' )
    doe.append( '\t'+'Custom Design,' )
    doe.append( '\t'+'{Add Response( Maximize, "Y", ., ., . ),' )
    npos = 0
    nfact = 0
    for pos in sorted(fact):
        name = fact[pos]['component']+str(pos)
        if len(fact[pos]['levels']) > 1:
            nfact += 1
            # Define as discrete no-empty factors (at least we need one!)
            if '-' not in fact[pos]['levels']:
                varType = 'Discrete Numeric'
                theLevels = [ x for x in range(1, len(fact[pos]['levels'])+1 ) ]
                doe.append( '\t'+'Add Factor( {varType}, {{{levels}}}, "{name}", 0),'.format(varType=varType,
                                                                                           levels=','.join(map(str,theLevels)),
                                                                                           name=name) )
            else:
                varType = 'Categorical'
                theLevels = [ '"L{}"'.format(x) for x in range(1, len(fact[pos]['levels'])+1 ) ]
                doe.append( '\t'+'Add Factor( {varType}, {{ {levels} }}, "{name}", 0),'.format(varType=varType,
                                                                                           levels=','.join(map(str,theLevels)),
                                                                                           name=name) )
        if fact[pos]['positional'] is not None:
            npos += 1
    if npos >  1:
        theLevels = ['"L{}"'.format(x) for x in range(1, npos*(npos-1)+1)]
        doe.append( '\t'+'Add Factor( {varType}, {{ {levels} }}, "{name}", 0),'.format(varType='Categorical',
                                                                                       levels=','.join(map(str,theLevels)),
                                                                                       name='pos') )
        nfact += 1
    if seed is None:
        seed = np.random.randint( 2**16 ) 
    doe.append( '\t'+'Set Random Seed( {} ),'.format( str(seed) ) )
    doe.append( '\t'+'Number of Starts( {} ),'.format( str(starts) ) )
    doe.append( '\t'+'Add Term( {1, 0} ),' )
    for i in range(1, nfact+1):
        doe.append('\t'+'Add Term( {{ {}, 1 }} ),'.format( str(i) ) )
    for i in range(1, nfact+1):
        for j in range(i+1, nfact+1):
            doe.append( '\t'+'Add Alias Term( {{ {}, 1 }}, {{ {}, 1 }} ),'.format( str(i), str(j) ))
    doe.append( '\t'+'Set Sample Size( {} ),'.format( str( int(size) ) ) )
    if makeTable:
        doe.append( '\t'+'Make Design,' )
        doe.append( '\t'+'Make Table' )
    else:
        doe.append( '\t'+'Make Design' )
    doe.append('\t'+'}')
    doe.append( ');')
    if makeTable:
        tabfile = re.sub("jsl$", "dat", outfile)
        doe.append( 'Current Data Table() << Save("{}", Text);'.format(tabfile) )
        doe.append( 'Quit();' )
    with open(outfile, 'w') as handler:
        handler.write( '\n'.join( doe ) )
                        

def addColumn(name, vartype, levels):
    par = []
    par.append( '"{}"'.format(name) )
    if vartype == "Discrete":
        par.append( 'Numeric')
        par.append( '"Continuous"' )
        par.append( 'Format( "Best", 12 )' )
        par.append( 'Set Property( "Design Role", Design Role( Discrete Numeric ) )' )
        par.append( 'Set Property( "Factor Changes", Easy )' )
#        for i in range(len(levels), totalLevels):
#            levels.append( '.' )
        values = '[{}]'.format( ','.join(map(str,levels)) )
        par.append( 'Set Values( {} )'.format(values) )
    elif vartype == "Nominal":
        par.append( 'Character' )
        par.append( '"Nominal"' )
        par.append( 'Set Property( "Design Role", Design Role( Categorical ) )' )
#        for i in range(len(levels), totalLevels):
#            levels.append( '.' )
        values = '{{"{}"}}'.format( '","'.join(levels) )
        par.append( 'Set Property( "Value Ordering", {} )'.format(values) )
        par.append( 'Set Property( "Factor Changes", Easy )' )
        par.append( 'Set Values( {} )'.format(values) )
    coldef = 'New Column( ' + ',\n\t'.join(par) + '\n)'
    return coldef
        

def makeTableScript(tableName, fact, outfile):
    """ Create a JMP script that generates the factor table """
    script = []
    npos = 0
    for pos in sorted(fact):
        name = fact[pos]['component']+str(pos)
        if len(fact[pos]['levels']) > 1:
            # Define as discrete no-empty factors (at least we need one!)
            if '-' not in fact[pos]['levels']:
                varType = 'Discrete'
                theLevels = [ x for x in range(1, len(fact[pos]['levels'])+1 ) ]
            else:
                varType = 'Nominal'
                theLevels = [ "L{}".format(x) for x in range(1, len(fact[pos]['levels'])+1 ) ]
            
            col = addColumn( name, varType, theLevels )
            script.append( col )
        if fact[pos]['positional'] is not None:
            npos += 1
    if npos >  1:
        theLevels = ["L{}".format(x) for x in range(1, npos*(npos-1)+1)]
        col = addColumn( 'pos', 'Nominal', theLevels )
        script.append( col )
    with open(outfile, 'w') as handler:
        handler.write( 'New Table( "{}",\n'.format( tableName ) )
        handler.write( 'Add Rows( {} ),\n'.format( len(script) ) )
        handler.write( 'New Table Variable( "Table Type", "DOE Factor Table" ),\n' )
        handler.write( ',\n'.join( script ) )
        handler.write( ')' )
    
def arguments():
    parser = argparse.ArgumentParser(description='Prepare JMP Files. Pablo Carbonell, SYNBIOCHEM, 2018')
    parser.add_argument('inputFile', 
                        help='Input file with specifications (excel format)')
    parser.add_argument('libSize', 
                        help='Library size')
    parser.add_argument('-O',  default=None,
                        help='Output folder (default: same as input)')
    parser.add_argument('-s', default=1,
                        help='Excel sheet number (default: 1)')
    parser.add_argument('-o',  default=None,
                        help='Output file (default: same as input file with jmp extension)')
    parser.add_argument('-r', action='store_true',
                        help='Overwrite if exists')
    parser.add_argument('-x', default=None,
                        help='Seed (default: random)')
    parser.add_argument('-n', default=100,
                        help='Number of starts (default: 100)')
    parser.add_argument('-l', default=None,
                        help='Log file') 
    parser.add_argument('-e', action='store_true',
                        help='Create an autoexecutable script')
    parser.add_argument('-t', action='store_true',
                        help='Make and save the table')
    return parser

if __name__ == '__main__':
    parser = arguments()
    arg = parser.parse_args()
    fact, seql, partinfo = read_excel( arg.inputFile )
    name = re.sub( '\.[^.]+$', '', os.path.basename(arg.inputFile) )
    outnametable = name+'_table.jsl'
    outname = name+'.jsl'
    logname = name+'.log'
    if arg.O is not None:
        if not os.path.exists( arg.O ):
            os.mkdir( arg.O )
        outdir = arg.O
    else:
        outdir = os.path.dirname(arg.inputFile)
    outfile = os.path.join( outdir, outnametable )
    if arg.o:
        outfile = arg.o
    if os.path.exists( outfile ) and not arg.r:
        raise Exception('File exists!')
    makeTableScript( name, fact, outfile )
    outfile = os.path.join( outdir, outname )
    if os.path.exists( outfile ) and not arg.r:
        raise Exception('File exists!')
    makeDoeScript( fact, outfile, size=arg.libSize, seed=arg.x, starts=arg.n, executable=arg.e, makeTable=arg.t )
    logfile = os.path.join( outdir, logname)
    with open(logfile, 'w') as handler:
        handler.write( ' '.join(['"{}"'.format(x) for x in sys.argv])+'\n' )
