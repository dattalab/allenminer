'''
Created on Jan 12, 2011

@author: andrewgiessel
'''

import matplotlib.pyplot as plt 
import numpy as np


if __name__ == '__main__':
    pass


def plotVectorField(X,Y):
    plt.figure()
    Q = plt.quiver(X,Y)
    qk = plt.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
                   fontproperties={'weight': 'bold'})
    l,r,b,t = plt.axis()
    dx, dy = r-l, t-b
    plt.axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])

    plt.title('Minimal arguments, no kwargs')
    plt.show()
    
def finiteDiffernceAxis(A,axis):
    X=np.empty_like(A)
    for coord in np.ndindex(A.shape):
        X[coord]=finiteDiffAtPoint(A,list(coord),axis)
    return X
            
def finiteDiffAtPoint(A,coord,axis):
    prevCoord=np.copy(coord)
    nextCoord=np.copy(coord)
    prevCoord[axis]-=1
    nextCoord[axis]+=1
    
    prevCoord=tuple(prevCoord)
    nextCoord=tuple(nextCoord)
    coord=tuple(coord)
    
    try:
        if prevCoord[axis]<0: raise IndexError  #prevent loop-around indexing
        prev=A[prevCoord]
    except:
        #use forward finite difference
        return A[nextCoord]-A[coord]
    try:
        next=A[nextCoord]
    except:
        #use backwards finite difference
        return A[coord]-A[prevCoord]
    
    #got both prev and next so do central difference
    return (next-prev)/2