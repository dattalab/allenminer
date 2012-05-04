'''
Created on Nov 29, 2010

@author: ag67
'''

import matplotlib as mplib
import matplotlib.pyplot as plt 
import numpy as np
import scipy as scipy
import scipy.stats as sstats
import os
import shelve
import re
import vectorfieldplot as vfp

import enthought.mayavi.mlab as mlab
from enthought.mayavi.api import *

if __name__ == "__main__":
    print "main"

def main():
    print "main"

#===============================================================================
#   Preprocessing, file reading and DB generating routines
#===============================================================================

def preprocessRoiListResults(filename):
    """Takes a roi_list output file, which has repeats of the following:
    headerline
    genebackground data (averaged across brain)
    gene data (for the roi of interest)
    
    This is read in and then split into two files, one for background 
    and one for the main genes.
    
    returns a tuple of string names of 2 new files if it worked.
    """
    
    gene_reader = open(filename, 'r')
    g_fh = open(("%s.genes" % filename), 'w')
    gb_fh = open(("%s.gene_backgrounds" % filename), 'w')
    for i, line in enumerate(gene_reader):
        if (np.mod(i, 3) == 0):     #skip header line
            continue
        elif (np.mod(i, 3) == 1):   #gene background line
            gb_fh.write(line)
        elif (np.mod(i, 3) == 2):   #gene line
            g_fh.write(line) 
    g_fh.close()
    gb_fh.close()
    gene_reader.close()
    return (''.join([filename, '.genes']), ''.join([filename, '.gene_backgrounds']))    

def killNaN(genes, gene_backgrounds):
    """Takes two numpy arrays.
    
    Deletes all rows where numpoints = 0 and returns the two cleaned arrays in a tuple
    """
    return (np.delete(genes,
                      [i for i, x in enumerate(genes['numpoints']) if x == 0.0],
                      0),
            np.delete(gene_backgrounds,
                      [i for i, x in enumerate(gene_backgrounds['numpoints']) if x == 0.0],
                      0))
        
def buildDetailedResultsDB(dbname, dirname, roipdbfilename):
    db = shelve.open(dbname)
    
    completeROI=readRoiPDB(roipdbfilename)
        
    for path, eachDir, files in os.walk(dirname):
        # parse filename into xpr name for key
        for eachFile in files:
            pat=re.compile('(.*)\.\d*\.roi.*')
            g=pat.match(eachFile).groups()
            key=g[0]
            print key
            try:
                result = readResult(path + '/' + eachFile)  #read in the np.array
                db[key] = detailedResult(key, result, completeROI)  # create a results object and put it in the DB
            except IOError:
                pass
    db.close()
    
def readRoiList(gene_filename, gene_backgrounds_filename):
    """Takes two strings which are names of gene and gene background data
    
    Returns two NumPy Record arrays of results-  genes and gene_backgrounds
    note: after stripping the NaNs out of these arrays
    """ 
    genes = np.genfromtxt(gene_filename, comments='#', delimiter='\t',
                     dtype=[('gene_name', '|S50'), ('slices', '|S50'), ('xpr_name', '|S50'), ('roi_name', '|S50'),
                            ('num_roidef_points', 'f8'), ('numpoints', 'f8'), ('expression_level', 'f8'),
                            ('nl_expression_level', 'f8'), ('num_expressors', 'f8'), ('nl_num_expressors', 'f8'),
                            ('specificity_expression_level', 'f8'), ('specificity_nl_expression_level', 'f8'), ('specificity_num_expressors', 'f8'),
                            ('specificity_nl_num_expressors', 'f8'), ('specificity_numpoints', 'f8')
                            ])
    gene_backgrounds = np.genfromtxt(gene_backgrounds_filename, comments='#', delimiter='\t',
                     dtype=[('gene_name', '|S50'), ('slices', '|S50'), ('xpr_name', '|S50'), ('roi_name', '|S50'),
                            ('num_roidef_points', 'f8'), ('numpoints', 'f8'), ('expression_level', 'f8'),
                            ('nl_expression_level', 'f8'), ('num_expressors', 'f8'), ('nl_num_expressors', 'f8'),
                            ('specificity_expression_level', 'f8'), ('specificity_nl_expression_level', 'f8'), ('specificity_num_expressors', 'f8'),
                            ('specificity_nl_num_expressors', 'f8'), ('specificity_numpoints', 'f8')
                            ])
    killDotXpr(genes)
    killDotXpr(gene_backgrounds)
    
    return killNaN(genes, gene_backgrounds)

def killDotXpr(g):
    for i, x in enumerate(g['xpr_name']):
        g[i]['xpr_name'] = x.split('.xpr')[0]

def readResult(filename):
    """Takes a filename of a detailed result hit.
    
    Returns a numpy rec array.
    
    Dtype= roi, x, y, z, expression_level, num_expressors, cell_diameter, grid_area
    """
    return np.genfromtxt(filename, comments='#', delimiter='\t',
             dtype=[('roi', '|S50'), ('x', 'i8'), ('y', 'i8'), ('z', 'i8'),
                    ('expression_level', 'f8'), ('num_expressors', 'i8'),
                    ('cell_diameter', 'f8'), ('grid_area', 'f8')])
    
def readRoiPDB(filename):
    """Takes a pdb filename of an ROI pdb.
    
    returns a tuple of np.arrays of X, Y and Z, for use with mlab.
    """  
    x = []
    y = []
    z = []
    pdb_reader = open(filename, 'r')
    for line in pdb_reader:
        row = line.split()
        try:
            x.append(float(row[6]))
            y.append(float(row[7]))
            z.append(float(row[8]))
        except IndexError:   
            pass
    pdb_reader.close()
    return (np.array(x), np.array(y), np.array(z))


#===============================================================================
#  plotting routines  (into object???)
#===============================================================================

def plotResult3d(result, scale_field='grid_area'):
    """Takes a detailed results numpy array
        
    plots the points, scaled by the num_expressors column
    """
    mlab.points3d(result.data['x'] - np.mean(result.data['x']),
                  result.data['y'] - np.mean(result.data['y']),
                  result.data['z'] - np.mean(result.data['z']),
                  result.data[scale_field])

def plotProjections(data):
    fig = plt.figure() 
    axy = fig.add_subplot(2, 2, 1)
    xyscat = plt.scatter(data['x'], data['y'])
        
    axz = fig.add_subplot(2, 2, 2)
    xzscat = plt.scatter(data['x'], data['z'])

    ayz = fig.add_subplot(2, 2, 3)
    yzscat = plt.scatter(data['y'], data['z'])

def plotROI3d(ROI):
    x = ROI['x']
    y = ROI['y']
    z = ROI['z']
    mlab.points3d(x - np.mean(x),
                  y - np.mean(y),
                  z - np.mean(z), opacity=.1)

def plotHit3d(results, n, engine=None, ROIcoord=None):
    if engine is not None:
        engine.scenes[0].children[1:2] = []
    if ROIcoord is not None:
        plotROI3d(ROIcoord)
    plotResult3d(results[n])
    print results[n].name
    
    # e=mlab.get_engine()
def plotPCs3d(v, s, PCAInputMatrix):
    #get extents
    longest_axis = max(PCAInputMatrix.max(axis=0) - PCAInputMatrix.min(axis=0))
    longest_axis /= 2
    scalar = s / s[0] * longest_axis
    for i, vector in enumerate(v):
        scaled_vector = vector * scalar[i]
        pc = zip(np.zeros(len(scaled_vector)), scaled_vector)
        mlab.plot3d(np.array(pc[0]), np.array(pc[1]), np.array(pc[2]), color=(1, 0, 0))

def plotImageArray(resultArray, firstCoord, secondCoord, numx, numy):
    fig=plt.figure()
    for i in range(numx*numy):
        a=fig.add_subplot(numx, numy, i+1)
        scat=plt.scatter(resultArray[i].data[firstCoord],resultArray[i].data[secondCoord])



#===============================================================================
# Stat/PCA routines
#===============================================================================

def centerOfMass(result, weightingField=None):
    """returns a tuple of the average of a result / ROI file.
    optionally weight by the indicated weighting field
    """
    
    if weightingField is not None:
        xmean=np.sum(result.data[weightingField]*result.data['x']) / np.sum(result.data[weightingField])
        ymean=np.sum(result.data[weightingField]*result.data['y']) / np.sum(result.data[weightingField])
        zmean=np.sum(result.data[weightingField]*result.data['z']) / np.sum(result.data[weightingField])
    else:
        xmean=np.mean(result.data['x'])
        ymean=np.mean(result.data['y'])
        zmean=np.mean(result.data['z'])
        try:
            xsem=sstats.sem(result.data['x'])
            ysem=sstats.sem(result.data['y'])
            zsem=sstats.sem(result.data['z'])
        except:
            xsem,ysem,zsem=(np.NaN,np.NaN,np.NaN)
    return xmean,ymean,zmean,xsem,ysem,zsem

def buildHist(resultList, dim):
    x,y,z=[],[],[]
    
    for eachResult in resultList:
        xmean,ymean,zmean,xsem,ysem,zsem=centerOfMass(eachResult)
        x.append(xmean)
        y.append(ymean)
        z.append(zmean)

    fig = plt.figure() 
    ax = fig.add_subplot(2, 2, 1)
    xhist = plt.hist(x,20)
        
    ay = fig.add_subplot(2, 2, 2)
    yhist = plt.hist(y,20)

    az = fig.add_subplot(2, 2, 3)
    zhist = plt.hist(z,20)


#===============================================================================
#   Util / PCA routines
#===============================================================================

def voxelNeigbors(coordinate):
    pass
#genrator to return nearest neighbors in 3d space
#    for i in [-1, 1]:
#        for j in [0,1,2]:
#            yield (coordinate[])

def coordListsFromPointTuples(tuples):
    x=[]
    y=[]
    z=[]
    for eachTuple in tuples:
        x.append(eachTuple[0])
        y.append(eachTuple[1])
        z.append(eachTuple[2])
    return (x,y,z)
    pass

def tuplesFromCoordLists(x,y,z):
    return [i for i in genFromCoordLists(x,y,z)]

def genFromCoordLists(x,y,z):
    for i in range(len(x)):
        yield int(x[i]), int(y[i]), int(z[i])
    
def pearsonsScore(roi, result1, result2):
    # for every point in the ROI, try to append the value in each
    # if not, then put a 0.0
        
    #blank flat arrays we're going to fill
    # will have one entry for each coord in the ROI
    # fill with the value from the complete orray, or 0.0 if it's null
    
    r1=[]
    r2=[]
    
    for coord in genFromCoordLists(roi['x'], roi['y'], roi['z']):
        for flat,result in zip([r1,r2],[result1,result2]):
            try:
                flat.append(result[coord]['grid_area'])
            except IndexError:
                flat.append(0.0)
    return sstats.pearsonr(r1, r2)
    
def buildPearsonCoMatrix(roi, resultList):
    N=len(resultList)

    m = scipy.zeros((N,N))
    
    for i, result1 in enumerate(resultList):
        for j, result2 in enumerate(resultList):
            p=pearsonsScore(roi, result1, result2)
            #print result1.name, result2.name, p
            m[i,j]=p[0]
    return m

    

def buildPCAInputMatrixFromResult(result, data='grid_area'):
    """ takes a result numpy array, optionally a value from the dtype list to put in the thing.
    
    returns a mean subtracted 4d m x n matrix.  m=4.  X, Y, Z and the data field
    n= number of expressing points.
    """
    xrow = result['x'] - np.mean(result['x'])
    yrow = result['y'] - np.mean(result['y'])
    zrow = result['z'] - np.mean(result['z'])
    datarow = result[data] - np.mean(result[data])
    
    return np.column_stack([xrow, yrow, zrow])

def calcPCA(result):
    result_data = buildPCAInputMatrixFromResult(result)
    u, s, v = scipy.linalg.svd(result_data)
    return result_data, u, s, v
    
def averageGradeint(resultArray):
    b=[]
    for eachResult in resultArray:
        a=eachResult.dataComplete()
        x=vfp.finiteDiff1DAxis(a, 0)
        y=vfp.finiteDiff1DAxis(a, 1)
        z=vfp.finiteDiff1DAxis(a, 2)
        b.append((eachResult.name, (x.mean(),y.mean(),x.mean())) )
    return b


#===============================================================================
# roi query result class defintion
#===============================================================================

class roiResult():
    def __init__(self, gene_filename, gene_background_filename, dbfilename, roipdbname):
        self.g, self.g_b = readRoiList(gene_filename, gene_background_filename)

        self.dbfilename = dbfilename
        self.roipdb = {}
        self.roipdb['x'], self.roipdb['y'], self.roipdb['z'] = readRoiPDB(roipdbname)
    
    def openResultsDB(self):    
        self.db = shelve.open(self.dbfilename)

    def closeResultsDB(self):
        self.db.close()
        
    def plotAll(self):
        """Plot all fields of a gene result array
        """ 
    
        fig = plt.figure() 
    
        for i in range(5, 14):
            ax = fig.add_subplot(9, 1, i - 4)
            ax.xlim = (0, 27000)
            l, = plt.plot(self.g[self.g.dtype.names[i]]) 
            t = ax.set_title(self.g.dtype.names[i]) 
    
    def sortgenes(self, columnToSort):
        """Takes a field name and sorts all genes based on this
    
        """
        self.g, self.g_b = np.sort(self.g, order=columnToSort), np.sort(self.g_b, order=columnToSort) 
        
    def plotNamed(self, columnToPlot):
        """Takes a string name of a field and plots it"""

        fig = plt.figure() 
        ax = fig.add_subplot(111)
        l, = plt.plot(self.g[columnToPlot])
              
    def plotROI(self):
        plotProjections(self.roipdb)
    
    def topHits(self, N):
        """Takes a number, N.
    
        Returns a list of detailed result objects.
        """
        self.openResultsDB()
        
        tophits = []
        L = len(self.g)
        for i in reversed(range(L - N, L)):
            name = self.g['xpr_name'][i]
            tophits.append(self.db[name])

        self.closeResultsDB()
        
        return tophits

    def roiExtent(self, coordinate, zeroMean=False):
        """Returns a tuple of the size of the ROI"""
        if not zeroMean:
            return (min(self.roipdb[coordinate]), max(self.roipdb[coordinate]))
        else:
            return (min(self.roipdb[coordinate])-np.mean(self.roipdb[coordinate]), 
                    np.max(self.roipdb[coordinate])-np.mean(self.roipdb[coordinate]))


class detailedResult():
    
    def __init__(self, name, data, roiPoints):
        self.name=name  #string of the name
        self.data=data  #numpy data
        
        self.min=[min(roiPoints[0]),
                  min(roiPoints[1]),
                  min(roiPoints[2])
                  ]
        
        self.max=[max(roiPoints[0]),
                  max(roiPoints[1]),
                  max(roiPoints[2])
                  ]

        self.dataDict = {}
        try:
            for row in self.data:
                self.dataDict[row['x'], row['y'], row['z']] = dict(zip(row.dtype.names[4:8], row.tolist()[4:8]))
        except TypeError:
            #we have a 1 line result.
            self.dataDict[int(self.data.tolist()[1]), int(self.data.tolist()[2]), int(self.data.tolist()[3])] = dict(zip(self.data.dtype.names[4:8], self.data.tolist()[4:8]))


        #dict of ndarrays of data, moved to the corner


    def dataComplete(self, dataField='grid_area'):       
        size=[]
        for coord in range(3):
            size.append(int(self.max[coord]-self.min[coord]+1))

        complete=np.zeros([size[0],size[1],size[2]])
        for coord, value in self.dataDict.iteritems():
            i=(coord[0]-self.min[0], coord[1]-self.min[1],coord[2]-self.min[2])
            complete[i]=value[dataField]
            
        return complete        
        
    def tupleInROI(self,coord):
        if (self.min[0]<=coord([0])<=self.max[0] and
            self.min[1]<=coord([1])<=self.max[1] and
            self.min[2]<=coord([2])<=self.max[2]):
            return True
        else:
            return False
        
    def __sub__(self, other):
        pass
    
    def __str__(self):
        return self.name+' detailed result object'
    
    def __getitem__(self, *args):
        try:
            return self.dataDict[args[0]]
        except:
            if tupleInRoi(args[0]):
                return {'cell_diameter': 0.0,
                        'expression_level': 0.0,
                        'grid_area': 0.0,
                        'num_expressors': 0.0}
            else:
                raise IndexError