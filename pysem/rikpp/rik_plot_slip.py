# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''
Main file to plot slip patches obtained with RIKsrf2

Example:
    
    python3 rik_plot_slip.py --nL 140 --nW 80 --sf /home/filippo/Data/Filippo/aeolus/ModeleVitessesEDF/slipdistribution.dat --hL 5.0 --hW 3.0 -L 7 -W 4
'''
#=======================================================================
# Required modules
#=======================================================================
import os
import sys
import argparse
import os.path as osp
from rik_pp_lib import *
#=======================================================================
# General informations
#=======================================================================
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

if __name__ == '__main__':
    opt = start_rik()
    globals().update(opt)
    # Input data
    NSR = nW*nL
    x, y, slip = readfileslip(NSR, sf)

    # Attention a l'ordre de W et L #
    plotslip (nW,nL,L,W,slip,True,fg,hL,hW)
