import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

class wavefunction():
    '''
    20240514: Indrasen Bhattacharya
    Class to explore complementarity in discrete signal processing
    Constructor initializes the operators in the chosen basis once
    Update can be used any number of times to change the wavefunction
    '''

    def __init__(self, initParams, funcParams, dtype_='complex64'):
        #Simple constructor
        #initialize using the position representation
        #N: the length of the array (number of sequential levels in the 'qudit')
        #M: (N-1)//2, N = 2*M+1, so that the total size is always odd
        #alpha: used in the inner product definition
        
        self.dtype = np.dtype(dtype_)

        self.N = initParams['N']
        self.M = (self.N-1)//2
        self.N = 2*self.M + 1
        self.alpha = initParams['alpha']
        self.initMethod = initParams['initMethod']
        self.indexArr = np.linspace(-self.M, self.M, self.N).astype('int')

        p_ = self.indexArr
        x_ = self.indexArr
        self.DFT = np.exp( -2j * np.pi * np.outer( x_, p_) / self.N ) / np.sqrt(self.N)

        self.Xoperator = self.initXop()
        self.Poperator = self.initPop()

        CC, AC = self.initCommutators()

        self.Commutator = CC
        self.Anticommutator = AC

        #this will update 
        #both the real space and momentum space descriptions
        self.update(funcParams)

        #math here is based on the indexArr being zero centered
        #manually construct the DFT matrix
        
        

    
    def update(self, funcParams):
        #update the wavefunction with a new set of function params
        #please note that all the functions are modified to have compact support
        #this leads to sinc convolution in the reciprocal space
        self.funcParams = {}
        self.initMethod =  funcParams['initMethod']
        self.funcParams['initMethod'] = self.initMethod #redundant

        if self.initMethod == 'box':
            
            self.funcParams['boxWidth'] = funcParams['boxWidth'] * self.N 
            self.funcParams['boxCenter'] = 0

            self.X = np.abs(self.indexArr - self.funcParams['boxCenter']) <= self.funcParams['boxWidth'] / 2 
            

        elif self.initMethod == 'gaussian':
            
            self.funcParams['sigma'] = funcParams['sigma'] * self.N
            self.funcParams['center'] = 0

            self.X = np.exp( -(self.indexArr - self.funcParams['center'])**2 / (2*(self.funcParams['sigma'])**2) ) 
            #self.X[ np.abs(self.indexArr) >= 3*self.funcParams['sigma'] ] = 0
            

        elif self.initMethod == 'lorentzian':
            #note the square root so that we can talk about normalization in the pdf sense
            
            self.funcParams['gamma'] = funcParams['gamma'] * self.N
            self.funcParams['center'] = 0

            self.X = np.sqrt( 1/(1 + ((self.indexArr - self.funcParams['center'])/self.funcParams['gamma'])**2) )
            #self.X[ np.abs(self.indexArr) >= 6*self.funcParams['gamma'] ] = 0

        elif self.initMethod == 'bump':
            #square root used here again

            self.funcParams['a'] = funcParams['a'] * self.N
            self.funcParams['center'] = 0

            self.X = np.sqrt( np.exp( -1 / (self.funcParams['a']**2 - (self.indexArr - self.funcParams['center'])**2) ) )
            self.X[ np.abs(self.indexArr - self.funcParams['center']) >= self.funcParams['a'] ] = 0

        elif self.initMethod == 'sinc':

            self.funcParams['a'] = funcParams['a'] * self.N
            self.funcParams['center'] = 0

            self.X = np.sin( np.pi * (self.indexArr - self.funcParams['center'] + np.finfo(np.float32).eps) / self.funcParams['a'] ) / (self.indexArr - self.funcParams['center'] + np.finfo(np.float32).eps) 
            #self.X[ np.abs( self.indexArr - self.funcParams['center'] ) >= 4*self.funcParams['a'] ] = 0

        elif self.initMethod == 'plug':
            #plug in the value from a certain input
            #will be useful when working with the eigenstates

            self.X = funcParams['psiX']
            
        
        self.X = np.asarray( self.X )
        self.X = self.X / np.power( self.innerProduct(self.X, self.X), 1/(2*self.alpha) )
        self.P = np.matmul( self.DFT, self.X )

        

    def innerProduct(self, vec1, vec2):
        return np.sum( (vec1**self.alpha) * np.conjugate(vec2**self.alpha) )
    
    def initXop(self):
        #Return the position operator in the position basis
        #A matrix is returned, expressed in the position basis

        return np.asarray( np.diag(self.indexArr), dtype=self.dtype)
    
    def initPop(self):
        #Return the momentum operator in the position basis
        #A matrix is returned, expressed in the position basis

        x = np.array( self.indexArr )
        diff = x[:, np.newaxis] - np.transpose( x[:, np.newaxis] )
        #r = np.exp( 2j * np.pi * diff / self.N )

        I = np.identity(self.N)
        A = np.ones_like(diff) - I

        #P = r**(-self.M) / (r-1) * (A - I)
        #P = r**(-self.M) * A
        Q = (-1)**(np.mod(diff,2)) * (2j) * np.sin( np.pi * diff / self.N )

        #diagonal entry doesn't matter since A makes it 0 anyway
        denom_ = np.reciprocal( Q + I )
        P = A * denom_


        return np.asarray( P, dtype=self.dtype )
    
    def initCommutators(self):
        #Return the commutator and anticommutator
        #Both are in the position basis

        x = np.array( self.indexArr )
        diff = x[:, np.newaxis] - np.transpose( x[:, np.newaxis] )

        x = np.array( self.indexArr, dtype=self.dtype )
        sum_ = x[:, np.newaxis] + np.transpose( x[:, np.newaxis] )

        I = np.identity(self.N, dtype=self.dtype)
        A = np.ones_like(diff) - I

        Q = (-1)**(np.mod(diff,2)) * (2j) * np.sin( np.pi * diff / self.N )

        #diagonal entry doesn't matter since A makes it 0 anyway
        denom_ = np.reciprocal( Q + I )

        CC = np.asarray( diff * A * denom_, dtype=self.dtype )
        AC = np.asarray( sum_ * A * denom_, dtype=self.dtype )

        return CC, AC
    
    def displ(self):
        #Display the probability densities in each of the spaces

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

        #density plot
        xProb = np.abs( (self.X ** self.alpha) * np.conjugate(self.X ** self.alpha) )
        pProb = np.abs( (self.P ** self.alpha) * np.conjugate(self.P ** self.alpha) )

        ax[0].set_title('PROBABILITY DISTRIBUTIONS')
        ax[0].plot(self.indexArr, xProb)
        ax[0].plot(self.indexArr, pProb)
        ax[0].legend(['POSITION SPACE', 'MOMENTUM SPACE'])
        ax[0].set_xlabel('INDEX')
        ax[0].set_ylabel('PROBABILITY DENSITY')

        #phase plot
        xPhase = np.angle( self.X )
        pPhase = np.angle( self.P )

        ax[1].set_title('PHASE DISTRIBUTION')
        ax[1].plot(self.indexArr, xPhase)
        ax[1].plot(self.indexArr, pPhase)
        ax[1].legend(['POSITION SPACE', 'MOMENTUM SPACE'])


    def returnUncertaintyProduct(self):
        #calculate using the position and momentum operators 
        #expressed in the position basis
        #remember that self.Poperator is expressed in the position basis
        #therefore, the expectation is taken with respect to the position basis wavefunction

        sigma_x = np.sqrt( np.abs(np.transpose( np.conjugate(self.X[:, np.newaxis]) ) @ self.Xoperator @ self.Xoperator @ (self.X[:, np.newaxis])) )
        sigma_p = np.sqrt( np.abs(np.transpose( np.conjugate(self.X[:, np.newaxis]) ) @ self.Poperator @ self.Poperator @ (self.X[:, np.newaxis])) )

        return np.squeeze( sigma_x * sigma_p ), np.squeeze(sigma_x), np.squeeze(sigma_p)

        

    def returnUncertaintyBound(self):
        #calculate using the commutator and anticommutator

        commutatorBound = np.abs(np.transpose( np.conjugate(self.X[:, np.newaxis]) ) @ self.Commutator @ (self.X[:, np.newaxis]))**2 / 4 
        anticommutatorBound = np.abs(np.transpose( np.conjugate(self.X[:, np.newaxis]) ) @ self.Anticommutator @ (self.X[:, np.newaxis]))**2 / 4 

        return np.squeeze( np.sqrt(commutatorBound + anticommutatorBound) )
    


def sweep(psiParams_, funcParams_, plotParams_, lowLim=0.001, upperLim=0.3):
    '''
    20240516: Indrasen Bhattacharya
    Function to sweep one parameter and generate the uncertainty plot
    '''

    psi = wavefunction(initParams=psiParams_, funcParams=funcParams_)

    commutatorBounds = []
    sigmaX = []
    sigmaP = []
    product = []

    sweep_arr = np.linspace(lowLim, upperLim, 100, endpoint=True)

    for sweepLoc in sweep_arr:
        funcParams_[plotParams_['qty']] = sweepLoc
        psi.update(funcParams_)
        
        commutatorBounds.append(psi.returnUncertaintyBound())
        
        prod, sigX, sigP = psi.returnUncertaintyProduct()

        product.append(prod)
        sigmaX.append(sigX)
        sigmaP.append(sigP)

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

    X = sweep_arr

    axs.semilogy(X, commutatorBounds)
    axs.semilogy(X, product)
    axs.semilogy(X, np.asarray(sigmaX))
    axs.semilogy(X, np.asarray(sigmaP))
    axs.legend(['BOUND', 'PRODUCT', 'STD(X)', 'STD(P)'])
    axs.grid(True)
    axs.set_xlabel(plotParams_['function'] + ' ' + plotParams_['qty'])
    axs.set_title('UNCERTAINTY PLOTS: ' + plotParams_['function'])
    