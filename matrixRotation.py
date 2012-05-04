#!/usr/bin/env python

'''
Created on Jan 30, 2011

@author: andrewgiessel
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import enthought.mayavi.mlab as mlab

def rot2d(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])

def rot3dX(theta):
    return np.array([[1,0,0], 
                     [0, np.cos(theta), -np.sin(theta)],
                     [0,np.sin(theta), np.cos(theta)]])
    
def rot3dY(theta):
    return np.array([[np.cos(theta), 0, np.sin(theta)], 
                     [0, 1, 0], 
                     [-np.sin(theta), 0, np.cos(theta)]])

def rot3dZ(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0], 
                     [np.sin(theta), np.cos(theta), 0], 
                     [0, 0,1]])


def make3dGradientKernel(size, low, high, thetaX, thetaY, thetaZ):
    # form data gradient

    krow=np.linspace(low, high, num=size)
    kplane=np.repeat(krow,size)
    k=np.repeat(kplane,size)

    # form coordinate meshes
    origCoords=np.mgrid[-(size/2):(size/2),-(size/2):(size/2),-(size/2):(size/2),]
    

    coordTuples=zip(origCoords[0].flatten(), origCoords[1].flatten(),origCoords[2].flatten())   # builds a list of coordinate tuples
    
    #these are the heart of the rotation.  uses nested maps (so read inside-out)
    # in the 3d form it is done by sucessive matrix mults.
    
    r=rot3dX(thetaX)
    rotCoordTuples = map(lambda y:(int(y[0]), int(y[1]), int(y[2])),               # turn into tuples of ints (for use in indexing)
                        map(np.around,                                 # round to nearest neighbor
                            map(lambda x:np.dot(x, r), coordTuples)))   # do the actual dot product of the coord and the rotation

    r=rot3dY(thetaY)
    rotCoordTuples = map(lambda y:(int(y[0]), int(y[1]), int(y[2])),               # turn into tuples of ints (for use in indexing)
                        map(np.around,                                 # round to nearest neighbor
                            map(lambda x:np.dot(x, r), rotCoordTuples)))   # do the actual dot product of the coord and the rotation

    r=rot3dZ(thetaZ)
    rotCoordTuples = map(lambda y:(int(y[0]), int(y[1]), int(y[2])),               # turn into tuples of ints (for use in indexing)
                        map(np.around,                                 # round to nearest neighbor
                            map(lambda x:np.dot(x, r), rotCoordTuples)))   # do the actual dot product of the coord and the rotation
    
# plot orig gradient and a translated rotated gradient (pre interp!) for comparison
    mlab.points3d(origCoords[0].flatten(), 
                  origCoords[1].flatten(), 
                  origCoords[2].flatten(),
                  k)
    mlab.points3d([c[0]+size*2 for c in rotCoordTuples],
                  [c[1] for c in rotCoordTuples],
                  [c[2] for c in rotCoordTuples],
                  k)

    rotatedKernel = buildKernelFromCoordTuplesAndData(rotCoordTuples,k) # build the translated X,Y matrix from the data and new coords
    rotatedKernel = interpolateData(rotatedKernel)  #fill the holes from rounding errors
    return rotatedKernel

def make2dGradientKernel(size, low, high, theta):
    # form data gradient
     
    krow = np.linspace(low, high, num=size)
    k = np.repeat(krow,size)
    
    # form coordinate meshes
    origCoords = np.mgrid[-(size/2):(size/2),-(size/2):(size/2)]

    #get the rotation matrix
    rrr = rot2d(theta)

    coordPairs = zip(origCoords[0].flatten(), origCoords[1].flatten())   # builds a list of coordinate tuples
    
    #this next line is the heart of the rotation.  uses nested maps (so read inside-out)
    
     
    rotCoordPairs = map(lambda y:(int(y[0]), int(y[1])),               # turn into tuples of ints (for use in indexing)
                map(np.around,                                 # round to nearest neighbor
                    rotCoordPairs = map(lambda x:np.dot(x, rrr), coordPairs)))   # do the actual dot product of the coord and the rotation

    1/0
    rotatedKernel = buildKernelFromCoordTuplesAndData(rotCoordPairs,k) # build the translated X,Y matrix from the data and new coords
    rotatedKernel = interpolateData(rotatedKernel)  #fill the holes from rounding errors
    return rotatedKernel
  
def buildKernelFromCoordTuplesAndData(coords, data):
    ndims=len(coords[0])
    mins=np.empty(ndims)
    maxes=np.empty(ndims)
    
    for c in range(ndims):
        coordinate=[yuyu[c] for yuyu in coords]
        maxes[c]=max(coordinate)
        mins[c]=min(coordinate)
    
    ranges=map(lambda x:x[0]*-1+x[1]+1, zip(mins,maxes) )
    
    fullarray=np.zeros(ranges)

    for c,datum in zip(coords,data):
        #have to offset the coords!
        offsets=map(lambda x:x[0]+x[1]*-1, zip(c,mins) )

        fullarray[tuple(offsets)]=datum
        
    return fullarray

#wrapper for 2d/3d interpolation
def interpolateData(rotatedKernel):
    if rotatedKernel.ndim==2:
        return interpolate2dData(rotatedKernel)
    if rotatedKernel.ndim==3:
        return interpolate3dData(rotatedKernel)


def interpolate2dData(rotatedKernel):
    #    check for holes!  it's not perfect interpolation but it's close-- and the larger the matrix the closer.  
    #    logic:  if point is zero (ie: never got assigned due to roundoff errors)
    #             and left and right neighbors are not, then set to average of left and right

    for i,row in enumerate(rotatedKernel):
        for j,point in enumerate(row):
                try:
                    if (j==0):                  # to get around pythons negative indexing "trick"
                        raise IndexError
                    prev=rotatedKernel[i,j-1]
                    next=rotatedKernel[i,j+1]
                except:
                    continue
                if (point==0 and prev!=0 and next!=0):   #are we a hole, inside core of the dataset?
                    rotatedKernel[i,j]=(prev+next)/2
    
    return rotatedKernel

def interpolate3dData(rotatedKernel):
    #    check for holes!  it's not perfect interpolation but it's close-- and the larger the matrix the closer.  
    #    logic:  if point is zero and left and right neighbors are not, then set to average of left and right

    for i,slice in enumerate(rotatedKernel):
        for j,row in enumerate(slice):
            for k,point in enumerate(row):
                try:
                    if (k==0):                  # to get around pythons negative indexing "trick"
                        raise IndexError
                    prev=rotatedKernel[i,j,k-1]
                    next=rotatedKernel[i,j,k+1]
                except:
                    continue
                if (point==0 and prev!=0 and next!=0):   #are we a hole, inside core of the dataset?
                    rotatedKernel[i,j,k]=(prev+next)/2
    
    return rotatedKernel

def makeAndPlotKernels(size,low,high,rotations):
    fig=plt.figure()
    
    for i, theta in enumerate(rotations):
        k=make2dGradientKernel(size, low, high, theta)
        ax = fig.add_subplot(1,len(rotations),i+1)
        plt.imshow(k)

#3d plot of a final rotated kernel
#unfort we have to go BACK to 3 lists of points + data list
def plotKernel3d(k):
    xs=[]
    ys=[]
    zs=[]
    datas=[]
    for x in range(k.shape[0]):
        for y in range(k.shape[1]):
            for z in range(k.shape[2]):
                xs.append(x)
                ys.append(y)
                zs.append(z)
                datas.append(k[x,y,z])
    mlab.points3d(xs,ys,zs,datas)



if __name__ == '__main__':
    print 'building kernels...'
    rotations=np.arange(0,16)/8*np.pi
    makeAndPlotKernels(100,0,10,rotations)
    plt.show()

