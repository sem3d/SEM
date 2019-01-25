# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to convert RIK output into SEM3D input
"""

#=======================================================================
# Required modules
#=======================================================================
from rik_pp_lib import *

#=======================================================================
# General informations
#=======================================================================
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, MSSMat UMR CNRS 8579 - CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

# Script to plot 3 subplots in a row

if __name__ == '__main__':
    fld = '/home/filippo/Data/Filippo/ares/workdir/RIK/napa_2014/results'
    # Input data
    Lgrid   = 150 
    Wgrid   = 100 
    LF = 15.
    WF = 10.
    NcX=12.5
    NcY=0.
    NSR = Wgrid*Lgrid

    filename = os.path.join(fld,'slipdistribution.dat')
    x, y, slip = readfileslip(NSR, filename)

    # Attention a l'ordre de W et L #
    plotslip (Wgrid,Lgrid,LF,WF,slip,False,'napa_slip.png',NcX,NcY)

    print '*** Done ***'
