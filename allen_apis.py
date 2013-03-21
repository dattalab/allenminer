import re
import array
import copy
import csv
import json
import struct
import sys
import urllib
import zipfile
import numpy as np


# lots of code here from the allen example.
# modified to return a np array instead of an array.array

# These are hard-coded paths to URLs for downloading expression volumes.
API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

STRUCTURE_GRAPH_ID = 1
REFERENCE_SPACE_ID = 10

STRUCTURES_URL = ("%s/OntologyNode/query.json?" +\
                      "criteria=[structure_graph_id$eq%d]") \
                      % (API_DATA_PATH, STRUCTURE_GRAPH_ID)

REFERENCE_SPACE_URL = ("%s/ReferenceSpace/query.json?criteria=[id$eq%d]" + \
                          "&include=well_known_files[path$li'*gridAnnotation.zip']" ) \
                          % (API_DATA_PATH, REFERENCE_SPACE_ID)

GRID_FMT = API_SERVER + "grid_data/download/%d"

# Make a query to the API via a URL.
def QueryAPI(url):
    start_row = 0
    num_rows = 2000
    total_rows = -1
    rows = []
    done = False

    # the ontology has to be downloaded in pages, since the API will not return
    # more than 2000 rows at once.
    while not done:
        pagedUrl = url + '&start_row=%d&num_rows=%d' % (start_row,num_rows)

        print pagedUrl
        source = urllib.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        
        if total_rows < 0:
            total_rows = int(response['total_rows'])

        start_row += len(response['msg'])

        if start_row >= total_rows:
            done = True

    return rows

# Download a grid file from the URL above by substituting in the data set id 
# argument.  Grid files are .zip files that will be downloaded to a 
# temporary location, where it can be unzipped into memory using the zipfile
# module.  The raw volume is converted into a flat array of floats.
def DownloadGridFile(dataSetId, dataset='energy'):
    assert dataset in ('energy', 'density', 'intensity', 'all'), 'DownloadGridFile: Unknown dataset \'%s\', must be one of \'energy\', \'density\', \'intensity\', \'all\'' % dataset
    
    if 'all' in dataset:
        url = (GRID_FMT + '?include=energy,intensity,density') % (dataSetId)
    elif 'energy' in dataset:
        url = (GRID_FMT  + '?include=energy') % (dataSetId)
    elif 'density' in dataset:
        url = (GRID_FMT + '?include=density') % (dataSetId) 
    elif 'intensity' in dataset:
        url = (GRID_FMT  + '?include=intensity') % (dataSetId)
       
    fh = urllib.urlretrieve(url)

    energy_header = None
    energy_data = None
    
    density_header = None
    density_data = None

    intensity_header = None
    intensity_data = None
        
    try:
        zf = zipfile.ZipFile(fh[0])

        if dataset is 'all' or dataset is 'energy':
            energy_header = zf.read('energy.mhd')
            energy_raw = zf.read('energy.raw')
            arr = array.array('f',energy_raw)
            energy_data = np.array(arr.tolist()).reshape(58, 41, 67)

        if dataset is 'all' or dataset is 'density':
            density_header = zf.read('density.mhd')
            density_raw = zf.read('density.raw')
            arr = array.array('f',density_raw)
            density_data = np.array(arr.tolist()).reshape(58, 41, 67)

        if dataset is 'all' or dataset is 'intensity':
            intensity_header = zf.read('intensity.mhd')
            intensity_raw = zf.read('intensity.raw')
            arr = array.array('f',intensity_raw)
            intensity_data = np.array(arr.tolist()).reshape(58, 41, 67)
    except:
        print "Bad zipfile for %s, returning null values!" % dataSetId
        return None
    
    
    data = {'energy_header':energy_header, 
            'energy_data':energy_data,
            'density_header':density_header,
            'density_data':density_data,
            'intensity_header':intensity_header,
            'intensity_data':intensity_data
            }     
        
    return data

# Download reference space meta information from the API.  Specifically, this 
# is looking for the download link to the zip file containing an annotation
# volume at the same resolution as the grid files.  Then, download the
# link, unzip the archive, and return the raw grid annotation volume as an array 
# of unsigned shorts (type `H`).
def DownloadAnnotationVolume():
    refspace = QueryAPI(REFERENCE_SPACE_URL)[0]
    reffile = refspace['well_known_files'][0]

    fh = urllib.urlretrieve(API_SERVER + reffile["download_link"])
    zf = zipfile.ZipFile(fh[0])

    raw = zf.read('gridAnnotation/gridAnnotation.raw')
    arr = array.array('H',raw)

    arr = np.array(arr.tolist()).reshape(58, 41, 67)
    
    return arr

# Download the ontology from the API.  This is a flat list of structures.
# We also download a list of parent-child relationships from the 
# StructureGraph model and use those to build a navigable tree hash.
def DownloadOntology():
    structures = QueryAPI(STRUCTURES_URL)

    # Build a hash from structure id to structure.
    structureHash = {}
    
    for i in xrange(len(structures)):
        s = structures[i]
        s['parent'] = None
        s['sum1'] = 0.0
        s['volume1'] = 0
        s['sum2'] = 0.0
        s['volume2'] = 0

        structureHash[s['structure_id']] = s

    # Make it a bit clearer who the parent of this structure is.
    for sid,s in structureHash.iteritems():
        if len(s['structure_id_path']) > 1:
            s['parent'] = structureHash[s['structure_id_path'][-2]]

    return structureHash
