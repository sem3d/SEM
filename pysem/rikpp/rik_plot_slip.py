# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''
Main file to plot slip patches obtained with RIKsrf2
'''
import os
import sys
import argparse
import os.path as osp
from rik_pp_lib import *
# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-L',type=float,default=15.,help="fault length")
    parser.add_argument('-W',type=float,default=15.,help="fault width")
    parser.add_argument('--nL',type=int,default=294,help="Number of grids along fault length")
    parser.add_argument('--nW',type=int,default=213,help="Number of grids along fault width")
    parser.add_argument('--hL',type=float,default=12.5,help="Hypocenter along-length coordinate")
    parser.add_argument('--hW',type=float,default=0.0,help="Hypocenter along-width coordinate")
    parser.add_argument('--sf',type=str,default="slipdistribution.dat",help="Filename of slip distributions")
    parser.add_argument('--fg',type=str,default="xxx.png",help="Figure file name")
    opt = parser.parse_args().__dict__
    globals().update(opt)
    # Input data
    NSR = nW*nL
    x, y, slip = readfileslip(NSR, sf)

    # Attention a l'ordre de W et L #
    plotslip (nW,nL,L,W,slip,False,fg,hL,hW)
