#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
learnFromData (c) University of Manchester 2019

learnFromData is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

Created on Thu Oct 31 14:37:54 2019

@author:  Pablo Carbonell, SYNBIOCHEM
@description: A stand-alone version of the learn part in mapPlasmids
"""
import os, re
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
import pylab as plt

option= 'InitialSelection' 
option = 'BestSelection' 
option = 'RegularSelection'
infile = '/mnt/SBC1/data/stress2/learn3/MAN13.xlsx'
infile2 = '/mnt/SBC1/data/stress2/learn/SBCDE00055_190814QQQMAN07_learn.csv'
infile3 = '/mnt/SBC1/data/stress2/learn/SBCDE00055_190814QQQMAN07_summary_0_constructs.csv'

def getFactors(df):
    factors = []
    for x in np.flip( df.columns ):
        if x == 'Design':
            break
        factors.append( x )
    factors.reverse()
    return factors

def getFormula(df,target,factors):
    formula = "Q('{}') ~".format( target )
    terms = []
    for f in factors:
        terms.append( "Q('{}')".format( f ) ) 
    formula += ' + '.join(terms)
    return formula

# Original values
df2 = pd.read_csv(infile2)
factors2 = getFactors(df2)
target2 = 'Target 1 Conc'
formula2 = getFormula(df2, target = target2, factors=factors2 )

# New values
df = pd.read_excel(infile)
df = df.loc[df.Design.notnull(),:]
if option == 'BestSelection':
    """ 
    In this case, 5 out of 6 top ranked are succesfully selected
    """
    target = 'Man_48'
    df = df.loc[df[target] > 10,:]
elif option == 'RegularSelection':
    target = 'Man_24'
       
# Keep only designs 55 and the 49 and 31 controls
df = df.loc[df.Design.str.contains(re.compile('(55)|(49)|(31)')),:]
#df = df.loc[df.Design != 'B5b',:]

# Comparison between experimental values

ex = dict( zip( df.Design, df[target] ))
ex2 = dict( zip( df2.Design, df2[target2] ))


common = list( set(ex) & set(ex2) )

x = [ex[z] for z in common]
y = [ex2[z] for z in common]

plt.scatter(y,x)
plt.xlabel('Original values')
plt.ylabel('New values')


# Build models
  
factors = getFactors(df)
formula = getFormula(df, target = target, factors=factors)

ols = smf.ols( formula=formula, data=df)
model = ols.fit()

ols2 = smf.ols( formula=formula2, data=df2)
model2 = ols2.fit()

# Comparison between fitted values

val = dict( zip( df.Design, model.fittedvalues ))
val2 = dict( zip( df2.Design, model2.fittedvalues ))

common = list( set(val) & set(val2) )

x = [val[z] for z in common]
y = [val2[z] for z in common]

plt.scatter(y,x)
plt.xlabel('Original values')
plt.ylabel('New values')

# Read constructs
df3 = pd.read_csv(infile3)
  
x = model.predict( df3 )
y = model2.predict( df3 )
plt.figure()

plt.scatter(y,x)
plt.xlabel('Predicted original values')
plt.ylabel('Predictd new values')

w1 = np.flip(np.argsort(x))
w2 = np.flip(np.argsort(y))
w = np.transpose(np.array([w1,w2,x[w1],y[w2]] ))
w2 = pd.DataFrame(w,columns=['New_rank','Old_rank', 'New_preds', 'Old_preds'])
print(w2.head(10))

## See how compare with predictions the new plasmids
# In[]

plt.figure()

if option == 'BestSelection':
    """ 
    In this case, 5 out of 6 top ranked are succesfully selected
    """
    df4 = pd.read_excel(infile)
    df4 = df4.loc[df4.Design.notnull(),:]
    df4 = df4.loc[df4.Design.str.contains('55'),:]

    target = 'Man_48'
    target2 = target
    df4 = df4.loc[df4[target] > 10,:]
elif option == 'RegularSelection':
    df4 = pd.read_excel(infile)
    df4 = df4.loc[df4.Design.notnull(),:]
    df4 = df4.loc[df4.Design.str.contains('55'),:]

    target = 'Man_24'
    target2 = target
    df4 = df4.loc[df4[target] > -1,:]

elif option == 'InitialSelection':
    df4 = pd.read_csv(infile2)
    df4 = df4.loc[df4.Design.notnull(),:]
    df4 = df4.loc[df4.Design.str.contains('55'),:]
    target = 'Target 1 Conc'
    target2 = 'Man_24'
else:
    assert False
    
factors4 = getFactors(df4)
formula4 = getFormula(df4, target = target, factors=factors4)
ols4 = smf.ols( formula=formula4, data=df4)
model4 = ols4.fit()



df5 = pd.read_excel(infile)
df5 = df5.loc[df5.Design.notnull(),:]
df5 = df5.loc[df5.Design.str.contains('56'),:]



df6 = pd.read_excel(infile)
df6 = df6.loc[df6.Design.notnull(),:]
df6 = df6.loc[df6.Design.str.contains('(49)|(31)'),:]



pred4 = model4.predict( df4 )

pred5 = model4.predict( df5 )
pred6 = model4.predict( df6 )

if option == 'InitialSelection':
    means = df4.groupby('Design').mean()[target]
    sds = df4.groupby('Design').std()[target]
    # To do: plot means + error bars
    plt.scatter(pred4, df4[target])

else:
    capsize = 3
    max_axis = 150
    fmt = 'D'
    means = df4[target]
    sdcol = target.split('_')
    sdcol = '_'.join( [sdcol[0],'sd',sdcol[1]] )
    #plt.scatter(pred4, df4[target])
#    fig = plt.figure()
#    ax = fig.add_axes([1,1,1,1])
#    ax.set_xlim([0,max_axis])
#    ax.set_ylim([0,max_axis])
    plt.xlim([0,max_axis])
    plt.ylim([0,max_axis])
    plt.errorbar(pred4,df4[target],yerr=df4[sdcol],fmt=fmt,capsize=capsize,color='royalblue')
    
    plt.errorbar(pred4,df4[target],yerr=df4[sdcol],fmt=fmt,capsize=capsize,color='royalblue')
    for i in range(0,len(df4.Name.index)):
        plt.text(pred4.iloc[i]+2,df4[target].iloc[i],s=df4.Name.iloc[i],color='royalblue')
#    plt.scatter(pred5, df5[target2],c='red')
    plt.errorbar(pred5,df5[target],yerr=df5[sdcol],fmt=fmt,capsize=capsize,color='green')
    for i in range(0,len(df5.Name.index)):
        plt.text(pred5.iloc[i]+2,df5[target].iloc[i],s=df5.Name.iloc[i],color='green')
#    plt.scatter(pred6, df6[target2],c='orange')
    means6 = df6.groupby('Design').mean()[target]
    sd6 = df6.groupby('Design').std()[target]
    plt.errorbar(pred6[0:2],means6,yerr=sd6,fmt=fmt,capsize=capsize,color='red')
    plt.text(pred6.iloc[0]+2,means6.iloc[0],s=df6.Name.iloc[0],color='red')
    plt.text(pred6.iloc[1]+2,means6.iloc[1],s=df6.Name.iloc[1],color='red')
#    plt.errorbar(pred6,df6[target],yerr=df6[sdcol],fmt='o')
    plt.xlabel('Predicted titers (mg/L)')
    plt.ylabel('Experimental titers (mg/L)')
#    plt.savefig('man_learn.png')
    plt.savefig('man_learn.svg')
#    plt.close(fig)
