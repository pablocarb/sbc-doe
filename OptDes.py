#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OptExp (c) University of Manchester 2018

OptExp is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>

Created on Tue Nov 27 16:01:46 2018

@author: pablo
"""
import numpy as np
import itertools

def Deff(X):
    # D-efficiency
    return (100.0/X.shape[0]) * ( np.linalg.det( np.dot( np.transpose( X ), X ) )**(1.0/X.shape[1]))
def Dopt(X):
    # D-optimality
    return np.linalg.det( np.dot( np.transpose( X ), X ) )
def SE(X):
    # Estimation efficiency
    return np.diag( np.linalg.inv( np.dot( np.transpose( X ), X ) ) )
def Contrib(X):
    cn = []
    for i in range(0, X.shape[0]):
        cn.append( Dopt( np.vstack( [X[:i,:], X[(i+1):,:]]  )  ) )
    return cn
def VarAdd(X,xj):
    # Variance of adding/removing one experiment
    return np.dot( np.dot( np.transpose(xj) , np.linalg.inv( np.dot( np.transpose( X ), X) ) ), xj )

def randExp( factors, n ):
    # Generate n random experiments
    V = None
    for levels in factors:
        vnew = np.random.randint(0, len(levels), n)
        if V is None:
            V = vnew
        else:
            V = np.vstack( [V, vnew] )
    return np.transpose( V )


#%%
def mapFactors( factors, M ):
    # Map a numerical factor into [-1,1] range, create dummy variables for categorical factors
    Mn = np.transpose( [np.ones( M.shape[0] )] )
    for i in np.arange( len(factors) ):
        v = factors[i]
        if type(v) == list:
            if len(set(v)) > 1:
                # Normalize between [-1,+1]
                Vn = (2*(M[:,i] - M[:,i].min())/(M[:,i].max()-M[:,i].min()) - 1)
                Vn = np.transpose( [Vn] )
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        else:
            if len(v) > 1:
                Vn = -np.ones( (M.shape[0],len(v)) )
                j = np.arange(M.shape[0])
                Vn[j,M[j,i]] = 1
                Vn = Vn[:,:-1]
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        if Mn is None:
            Mn = Vn
        else:
            Mn = np.hstack( [Mn, Vn])
    return Mn

#%%

def fullFactorial( factors ):
    # Here we generate a full factorial but this is not possible for large designs
    # Replace by random sampling, ideally some descent algorithm (TBD)
    # For categorical variables, we randomize the levels that are then mapped into the dummy variables
    val = []
    # Add a constant for the intercept
    val.append( [1] )
    for v in factors:
        val.append( np.arange(len(v)) )
    ffact = []
    for m in itertools.product(*val):
        ffact.append(m)
    ffact = np.array(ffact)
    return ffact

#%%


def DetMax( factors, n, m, w=1000, th=99.5, k=1 ):
    # n: Number of runs
    # m: sampled design space per iteration
    # w: maximum iterations
    # k: number of removal/additions per iteration
    # th: stop criterium for D-efficiency (<= 100.0)
    # Here I have implemented a simple DETMAX algorithm. At each iteration: 
    #  - remove the design with the lowest variance
    # - add the design with the highest variance
    # Many more efficent variants exist (kl-exchange, genetic algorithms, etc..)
    # See Mandal et al, Algorithmic searches for optimal design
    
    # Initial design: could we do something better than a purely random start?   
    M = randExp( factors, n )
    # Map factors into [-1,1] and convert chategorical
    X = mapFactors( factors, M )

    # D-Efficiency of the initial design
    J = Deff(X)
    print(J)
    w = 0
    while ((J<100.0) and (w < 1000)):
        # X1 is the design space sample in the iteration. 
        # Here we loop through the full factorial, which is computationally expensive
        # First thing to fdo is to change it to random generation of a subset of the full library
        # It would be better to move across some surface like gradient descent...
        M1 = randExp( factors, m )
        X1 = mapFactors( factors, M1 ) 
        # Calcualte delta of removing an experiment
        sub = []
        for i in np.arange(X.shape[0]):
            sub.append( VarAdd(X, X[i,:]) )
        w += 1
        # Remove the experiments with the lowest delta
        Xsub = None
        dList = np.argsort( sub )[0:(k-1)]
        for i in np.arange(X.shape[0]):
            if i in dList:
                continue
            else:
                if Xsub is None:
                    Xsub = X[i,:]
                else:
                    Xsub = np.vstack( [Xsub, X[i,:]] )
        # Calculate the delta of adding an experiment from the sample
        add = []
        for j in np.arange(X1.shape[0]):
            add.append( VarAdd( Xsub, X1[j,:] ) )
        # Add the experiments with the highest delta
        aList = np.flip( np.argsort( add ) )[0:(k+1)]
        Xn = Xsub
        for j in aList:
            Xn = np.vstack( [Xn, X1[j,:] ] )
        # Make the update id the resulting design increases the objective
        #print(w,J,i,j, Dopt(X), Dopt(Xsub), Dopt(Xn))
        if Dopt(Xn) > Dopt(X):
            X = Xn
            J = Deff(X)
        elif Dopt(Xn) == Dopt(X):
            break
    print(w,J)
    return X
    

n = 46 # Number of runs
m = 100 # Sampled design space per iteration

factors = [ [0,1,2,3,4], 
           [123,53,345],#[0,1], [0,1], [0,1]]
   {'Red', 'Green', 'Blue'}, 
   {'prom1', 'prom2', 'prom3', 'prom4'} ]

# Initial design: could we do something better than a purely random start?
    
M = randExp( factors, n )
X = mapFactors( factors, M )

# D-Efficiency of the initial design
# Here I have implemented a simple DETMAX algorithm. At each iteration: 
#  - remove the design with the lowest variance
#  - add the design with the highest variance
# Many more efficent variants exist (kl-exchange, etc..)
J = Deff(X)
print(J)
w = 0
while ((J<100.0) and (w < 1000)):
    # X1 is the design space sample in the iteration. 
    # Here we loop through the full factorial, which is computationally expensive
    # First thing to fdo is to change it to random generation of a subset of the full library
    # It would be better to move across some surface like gradient descent...
    M1 = randExp( factors, m )
    X1 = mapFactors( factors, M1 ) 
    sub = []
    for i in np.arange(X.shape[0]):
        sub.append( VarAdd(X, X[i,:]) )
    w += 1
    Xsub = None
    dList = np.argsort( sub )[0:1]
    for i in np.arange(X.shape[0]):
        if i in dList:
            continue
        else:
            if Xsub is None:
                Xsub = X[i,:]
            else:
                Xsub = np.vstack( [Xsub, X[i,:]] )
    add = []
    for j in np.arange(X1.shape[0]):
        add.append( VarAdd( Xsub, X1[j,:] ) ) 
    aList = np.flip( np.argsort( add ) )[0:1]
    Xn = Xsub
    for j in aList:
        Xn = np.vstack( [Xn, X1[j,:] ] )
#    if w % 10 == 0:
    print(w,J,i,j, Dopt(X), Dopt(Xsub), Dopt(Xn))
    if Dopt(Xn) > Dopt(X):
        X = Xn
        J = Deff(X)
    elif Dopt(Xn) == Dopt(X):
        break
print(w,J)

# Define two type of factors:
# - Numeric (discrete?) => ordered list
# - Categorical => set


