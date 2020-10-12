# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
XDMF helper
"""
#===============================================================================
# Required modules
#===============================================================================
# General modules
from __future__ import division
import warnings
warnings.filterwarnings("ignore")
import os.path as osp
import argparse
import numpy as np
import h5py
#===============================================================================
# General informations
#===============================================================================
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

def write_header(fid):
    fid.write('''<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n'''+
    '''<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n'''+
    '''<Domain>\n''')
    fid.write('''<Grid Name="{}" GridType="Collection" CollectionType="Temporal">\n'''.format('fault'))
    return

def write_close(fid):
    fid.write('''</Grid>\n</Domain>\n</Xdmf>\n''')
    return

def write_geometry(fid,h5f,tfm,nel=1,nnd=1,tag='Box',gtp='XYZ'):
    
    fid.write('''<Grid Name="{}">\n'''.format(tag))
    fid.write('''<Time Value="{:<20.10f}"/>\n'''.format(tfm))
    fid.write('''<Topology TopologyType="Triangle" NumberOfElements="{:>d}">\n'''.format(nel))
    fid.write('''<DataItem Format="HDF" NumberType="Int" Dimensions="{:>d}  3">\n'''.format(nel))
    fid.write('''{}:/Elements\n'''.format(h5f))
    fid.write('''</DataItem>\n''')
    fid.write('''</Topology>\n''')
    fid.write('''<Geometry GeometryType="{}">\n'''.format(gtp))
    fid.write('''<DataItem Format="HDF" NumberType="Float" Precision="4"'''+
              ''' Dimensions=" {:>15} 3">\n'''.format(nnd)+
              '''{}:/Nodes\n'''.format(h5f)+
              '''</DataItem>\n'''+
              '''</Geometry>\n''')
    return

def write_attributes(fid,h5f,attrs):
    for a in attrs:
        for k,v in a.items():
            fid.write('''<Attribute Name="{}" Center="{}" AttributeType="{}">\n'''.format(v[0],v[1],v[2]))
            if len(v[3])>1:
                fid.write('''<DataItem Format="HDF" NumberType="Float" Precision="4"'''+\
                          ''' Dimensions=" {} {}">\n'''.format(*v[3]))
            else:
                fid.write('''<DataItem Format="HDF" NumberType="Float" Precision="4"'''+\
                          ''' Dimensions=" {}">\n'''.format(*v[3]))
            fid.write('''{}:/{}</DataItem>\n</Attribute>\n'''.format(h5f,k))
    return

