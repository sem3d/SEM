# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''
Main file to plot slip patches obtained with RIKsrf2

Example:

    Plot slip contour over 140 x 80 fault grid, with hypocenter located at (x,y=(5 km, 3 km) on the fault 2D plane 
    ^ W 
    |
    |____> L

    python3 rik_plot_slip.py @wkd ./ @nL 140 @nW 80 @sf ./slipdistribution.dat @hL 3.66 @hW 3.15 @L 7.0 @W 4.0 @tag teil

'''
#=======================================================================
# Required modules
#=======================================================================
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
    fg=opj(wkd,tag)
    NSR = nW*nL
    x, y, slip = readfileslip(NSR, sf)
    # Attention a l'ordre de W et L #
    plotslip (nW,nL,L,W,slip,True,fg,hL,hW)
