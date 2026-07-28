# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to create shapes for new materials
"""
#===============================================================================
# Required modules
#===============================================================================
# General modules
import h5py
import numpy as np
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

class sphere(object):
    def __init__(self,xc=0.,yc=0.,zc=0.,Rd=1.):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.Rd = Rd
        self.volume()
    def __call__(self,xc=0.,yc=0.,zc=0.,Rd=1.):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.Rd = Rd
        self.vol = 4./3.*np.pi*Rd**3;
    def iswithin(self,x,y,z):
        ixyz = ((x-self.xc)**2+(y-self.yc)**2+(z-self.zc)**2 <= self.Rd**2)
        return ixyz
    def volume(self):
        self.vol = 4./3.*np.pi*self.Rd**3;
        return


class ellipsoid(object):
    def __init__(self,xc=0.,yc=0.,zc=0.,Ad=4.5,Bd=6.0,Cd=3.0):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.Ad = Ad
        self.Bd = Bd
        self.Cd = Cd
        self.volume()
    def __call__(self,xc=0.,yc=0.,zc=0.,Ad=4.5,Bd=6.0,Cd=3.0):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.Ad = Ad
        self.Bd = Bd
        self.Cd = Cd
        self.vol = 4./3.*np.pi*self.Ad*self.Bd*self.Cd;
    def iswithin(self,x,y,z):
        ixyz = ((x-self.xc)**2/(self.Ad**2)+\
                (y-self.yc)**2/(self.Bd**2)+\
                (z-self.zc)**2/(self.Cd**2) <= 1.0)
        return ixyz
    def volume(self):
        self.vol = 4./3.*np.pi*self.Ad*self.Bd*self.Cd;
        return
