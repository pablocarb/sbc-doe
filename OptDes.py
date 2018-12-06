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
import pandas as pd
import itertools

def Deff(X):
    # D-efficiency
    return (100.0/X.shape[0]) * ( np.linalg.det( np.dot( np.transpose( X ), X ) )**(1.0/X.shape[1]))

def Deff2(M, factors):
    X = mapFactors2(M, factors)
    return (100.0/X.shape[0]) * ( np.linalg.det( np.dot( np.transpose( X ), X ) )**(1.0/X.shape[1]))

def Dopt(X):
    # D-optimality
    return np.linalg.det( np.dot( np.transpose( X ), X ) )

def Dopt2(M, factors):
    # D-optimality
    X = mapFactors2(M, factors)
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
    
def grid(n, weighted=True):
    """ Provide normalized vectors of n-1 dummy variables 
    Useful for computing the model matrix (X) to 
    use pseudo-orthogonal terms in the n-1 hypercube.
    (Experimental)
    In JMP, grid(3) is multiplied by sqrt(2), grid(4) by
    sqrt(3), which brings back the weight of the number of 
    factors
    """
    from sklearn.preprocessing import normalize
    from sklearn.decomposition import PCA
    
    base = np.eye(n)*2 - 1
    pc = PCA(n-1, whiten=True, random_state=0)
    bt = pc.fit_transform(base)
    W = normalize(bt)
    if weighted:
        W = W*np.sqrt(n-1)
    return W

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

def mapFactors2( M, factors ):
    # Map a numerical factor into [-1,1] range, 
    # create orthogonal coordinates for dummy variables for categorical factors
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
                # Use grid 
                j = np.arange(M.shape[0])
                Vn = gridList[len(v)][M[j,i],:]
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        if Mn is None:
            Mn = Vn
        else:
            Mn = np.hstack( [Mn, Vn])
    return Mn


def MapExp( E ):
    """ Read a design, transform into X matrix """
    # Define factors in the same way as for the library
    factors = [ set(np.unique(E[:,i])) for i in np.arange(E.shape[1])]
    EE = np.transpose( np.array([ list(np.unique(E[:,i], return_inverse=True)[1]) for i in np.arange(E.shape[1])] ) )
    M = mapFactors( factors, EE )
    return M, factors, EE

def MapDesign(factors, X):
    """ Map back from X to the factors """
    M = []
    for i in np.arange(X.shape[0]):
        row = []
        # skip intercept
        j = 1
        for fa in factors:
            levels = sorted(fa)
            # If none is set
            level = levels[-1]
            for l in levels[0:-1]:
                if X[i,j] == 1:
                    level = l
                j += 1
            row.append(level)
        M.append( row )
    return np.array( M )
    
def MapDesign2(factors, M):
    """ Map back from M to the factors """
    N = []
    for i in np.arange(M.shape[0]):
        row = []
        for j in np.arange(M.shape[1]):
            levels = sorted(factors[j])
            row.append( levels[ M[i,j] ])
        N.append( row )
    return np.array( N )


def JMPExample():
    """ This is a JMP example: """
    # Design Evaluation
    # Design Diagnostics
    # D Optimal Design	
    # D Efficiency	 87.98414
    # G Efficiency	 64.62616
    # A Efficiency	 76.00696
    # Average Variance of Prediction 1.229865
    # Design Creation Time (seconds) 11
    # Read design
    E = pd.read_excel('/mnt/SBC1/data/OptimalDesign/data/CD2.xlsx')
    D = np.array(E.iloc[:,0:8])
    # Map into a binary matrix 
    DD, fac, EE = MapExp(D)
    # Compute D-efficiency (wrong use of categorical factors)
    print( Deff( DD ) )
    # 38.66
    # Compute D-efficiency (correct)
    print( Deff2( EE, factors) )
    # D Efficiency	 87.98414
    return fac, DD, EE

    
    
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


def DetMax( factors, n, m, it=1000, th=99.5, k=1 ):
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
    # Map factors into [-1,1] and convert categorical
    X = mapFactors( factors, M )
    # D-Efficiency of the initial design
    J = Deff(X)
    print(J)
    w = 0
    while ((J<100.0) and (w < it)):
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
        dList = np.argsort( sub )[0:(k+1)]
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
        aList = np.flip( np.argsort( add ), axis=0 )[0:(k+1)]
        Xn = Xsub
        for j in aList:
            Xn = np.vstack( [Xn, X1[j,:] ] )
        # Make the update if the resulting design increases the objective
        if w % 10 == 0:
            print(w,J,i,j, Dopt(X), Dopt(Xsub), Dopt(Xn))
        if Dopt(Xn) > Dopt(X):
            X = Xn
            J = Deff(X)
        elif Dopt(Xn) == Dopt(X):
            break
    print(w,J)
    return X

def DetMax2( factors, n, m, it=1000, th=99.5, k=1, verb=False ):
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
    X = mapFactors2( M, factors ) 
    # D-Efficiency of the initial design
    J = Deff2(M, factors)
    if verb:
        print(J)
    w = 0
    while ((J<100.0) and (w < it)):
        # X1 is the design space sample in the iteration. 
        # Here we loop through the full factorial, which is computationally expensive
        # First thing to fdo is to change it to random generation of a subset of the full library
        # It would be better to move across some surface like gradient descent...
        M1 = randExp( factors, m )
        X1 = mapFactors2( M1, factors ) 
        # Calcualte delta of removing an experiment
        sub = []
        for i in np.arange(X.shape[0]):
            sub.append( VarAdd(X, X[i,:]) )
        w += 1
        # Remove the experiments with the lowest delta
        Xsub = None
        dList = np.argsort( sub )[0:(k+1)]
        for i in np.arange(X.shape[0]):
            if i in dList:
                continue
            else:
                if Xsub is None:
                    Xsub = X[i,:]
                    Msub = M[i,:]
                else:
                    Xsub = np.vstack( [Xsub, X[i,:]] )
                    Msub = np.vstack( [Msub, M[i,:]] )
        # Calculate the delta of adding an experiment from the sample
        add = []
        for j in np.arange(X1.shape[0]):
            add.append( VarAdd( Xsub, X1[j,:] ) )
        # Add the experiments with the highest delta
        aList = np.flip( np.argsort( add ), axis=0 )[0:(k+1)]
        Xn = Xsub
        Mn = Msub
        for j in aList:
            Xn = np.vstack( [Xn, X1[j,:] ] )
            Mn = np.vstack( [Mn, M1[j,:] ] )
        # Make the update if the resulting design increases the objective
        if w % 10 == 0 and verb:
            print(w,J,i,j, Deff2(M,factors))
        if Dopt(Xn) > Dopt(X):
            X = Xn
            M = Mn
            J = Deff2(M,factors)
        elif Dopt(Xn) == Dopt(X):
            break
    if verb:
        print(w,J)
    return M


def blending(A,B):
    blend = np.random.randint(2, size=A.shape[0])
    ib = np.argwhere( blend == 1 )
    offspring = A.copy()
    offspring[ib] = B[ib]
    return offspring

def crossover(A,B):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    x1 = A[i,np.arange(0,j)]
    x2 = A[i,np.arange(j,A.shape[1])]
    y1 = B[i,np.arange(0,j)]
    y2 = B[i,np.arange(j,A.shape[1])]
    A[i] = np.append(x1,y2)
    B[i] = np.append(y1,x2)

def mutation(A, th=2.0):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    epsilon = np.random.normal()
    if epsilon > th:
        if A[i,j] > 0:
            A[i,j] = -1
        else:
            A[i,j] = 1
        

def reproduction(A,B):
    A = A.copy()
    B = B.copy()
    for i in np.arange(100):
        if np.random.randint(100) > 95:
            mutation(A)
        if np.random.randint(100) > 95:
            mutation(B)
    # Crossover for multilevel numerical factors won't work
    for i in np.arange(10):
       if np.random.randint(100) > 80:
            crossover(A,B)
    C = blending(A,B)
    return C
    

    
def GenAlg( factors, n, m, nPop=100, it=10, th=99.5, k=1, func=Dopt):
    # Using genetic algorithms
    # Based on Using a Genetic Algorithm to Generate Doptimal Designs for Mixture Experiments
    # To be revised: mutations need to be done in the multilevel design space
    # 1. Generate an initial population of random designs
    # Take the best ones
    population = []
    for i in np.arange(nPop*20):
        M = randExp( factors, n )
        X = mapFactors( factors, M )
        population.append( X )
    population = np.array(population)
    eff = []
    for j in np.arange(nPop):
        eff = np.append( eff, func(population[j]) )
    effi = np.flip( np.argsort(eff), axis=0)
    population = population[effi[0:nPop]]
    # 2. Calculate score for each member of the population
    # and select elite chromosome
    for i in np.arange(it+1):
        eff = np.array( [] )
        for j in np.arange(nPop):
            eff = np.append( eff, func(population[j]) )
        effi = np.flip( np.argsort(eff), axis=0)
        elite = effi[0]
        print(i,Deff(population[elite]))
        if i == it:
            return population[elite]
        population = population[effi[0:nPop]]
        eff = eff[effi[0:nPop]]
#        population = np
        # 3. Random pairing
#        w = np.arange(1, len(eff)-1)
#        np.random.shuffle( eff[w] )
        pairs = []
  #      for i in np.arange(len(w), step=2):
        for i in np.arange(50):
            a = np.random.randint(50)
            b = np.random.randint(20)
            if a != b:
                pairs.append( (population[a], population[b]) )
        for p in pairs:
            offspring = reproduction(p[0],p[1])
            population = np.insert(population,-1,offspring, axis=0)
    
def blending2(A,B):
    blend = np.random.randint(2, size=A.shape[0])
    ib = np.argwhere( blend == 1 )
    offspring = A.copy()
    offspring[ib] = B[ib]
    return offspring

def crossover2(A,B):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    x1 = A[i,np.arange(0,j)]
    x2 = A[i,np.arange(j,A.shape[1])]
    y1 = B[i,np.arange(0,j)]
    y2 = B[i,np.arange(j,A.shape[1])]
    if np.random.randint(2)>0:
        A[i] = np.append(x1,y2)
        B[i] = np.append(y1,x2)
    else:
        A[i] = np.append(y1,x2)
        B[i] = np.append(x1,y2)

def mutation2(A, th=1.0):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    epsilon = np.random.normal()
    if epsilon > th:
        levels = np.unique(A[:,j])
        A[i,j] = levels[ np.random.randint(len(levels)) ]
        
def reproduction2(A,B):
    A = A.copy()
    B = B.copy()
    for i in np.arange(100):
        mutation2(A)
        mutation2(B)
    # Crossover for multilevel numerical factors won't work
    for i in np.arange(5):
       if np.random.randint(100) > 50:
            crossover2(A,B)
    C = blending2(A,B)
    return C
  
    
def GenAlg2( factors, n, nPop=100, it=10, th=99.5 ):
    # Using genetic algorithms
    # Based on Using a Genetic Algorithm to Generate Doptimal Designs for Mixture Experiments
    # To be revised: mutations need to be done in the multilevel design space
    # 1. Generate an initial population of random designs
    # Take the best ones
    population = []
    for i in np.arange(nPop):
#        M = randExp( factors, n )
        M = DetMax2(factors, n, m=100, it=10, k=2)
        population.append( M )
    population = np.array(population)
    eff = []
    for j in np.arange(nPop):
        eff = np.append( eff, Deff2(population[j], factors) )
    eff[ np.isnan(eff) ] = 0
    effi = np.flip( np.argsort(eff), axis=0)
    population = population[effi[0:nPop]]
    # 2. Calculate score for each member of the population
    # and select elite chromosome
    for i in np.arange(it+1):
        eff = np.array( [] )
        for j in np.arange(len(population)):
            eff = np.append( eff, Deff2(population[j], factors) )
        eff[ np.isnan(eff) ] = 0
        effi = np.flip( np.argsort(eff), axis=0)
        elite = effi[0]
        print(i,Deff2(population[elite],factors),Deff2(population[effi[1]],factors))
        if i == it:
            return population[elite]
        population = population[effi[0:nPop]]
        eff = eff[effi[0:nPop]]
#        population = np
        # 3. Random pairing
#        w = np.arange(1, len(eff)-1)
#        np.random.shuffle( eff[w] )
        pairs = []
  #      for i in np.arange(len(w), step=2):
        for a in np.arange(10):
            for b in np.arange(10):
                if a != b:
                    for k in np.arange(10):
                        pairs.append( (population[a], population[b]) )
        for p in pairs:
            offspring = reproduction2(p[0],p[1])
            population = np.insert(population,-1,offspring, axis=0)
#%%    

# Precompute the hypercube grids
gridList = {}
for i in np.arange(2,20):
    gridList[i] = grid(i)
    
n = 46 # Number of runs
m = 100 # Sampled design space per iteration

factors = [ [0,1,2,3,4], 
           [123,53,345],#[0,1], [0,1], [0,1]]
   {'Red', 'Green', 'Blue'}, 
   {'prom1', 'prom2', 'prom3', 'prom4'} ]

factors = [
   {'Red', 'Green', 'Blue'}, 
   {'prom1', 'prom2', 'prom3', 'prom4','prom1', 'prom2', 'prom3', 'prom4'},
   {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
   {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
   {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
   {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
   {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
   {'prom1', 'prom2', 'prom3', 'prom4'} ]

factors,DD, EE = JMPExample()

#X = DetMax(factors, n, m, it=1000, k=2)


#M = DetMax2(factors, n, m, it=100, k=2)

best = 0
while best < 100:
#    M = DetMax2(factors, n, m=100, it=10, k=2)
    M = GenAlg2(factors, n=46, it=10)
    de = Deff2(M, factors)
    
    if de > best:
        best= de
        print(best)


#X = GenAlg(factors, n, m=20, it=100, func=Dopt)

#X = GenAlg2(factors, n=46, it=100)

