# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Set of subroutines to post-process RIK output
"""
# Required modules
import utm
import os 
from os.path import join as opj
import argparse
import sys
import numpy as np
import math
import shutil
import h5py
from scipy import signal
from scipy.interpolate import interp1d
import pylab as p
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib import ticker
from matplotlib.colors import LogNorm
from matplotlib.image import NonUniformImage
from matplotlib.ticker import LogFormatter, FormatStrFormatter
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm,rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter
rc('text', usetex=True)
import seaborn as sns
from   collections import OrderedDict 

# General informations
__author__ = "Filippo Gatti, Elif Oral"
__copyright__ = "Copyright 2023, LMPS UMR 9026, CentraleSupélec"
__credits__ = ["Filippo Gatti", "Elif Oral"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"


class ClassPropertyDescriptor(object):

    def __init__(self, fget, fset=None):
        self.fget = fget
        self.fset = fset

    def __get__(self, obj, klass=None):
        if klass is None:
            klass = type(obj)
        return self.fget.__get__(obj, klass)()

    def __set__(self, obj, value):
        if not self.fset:
            raise AttributeError("can't set attribute")
        type_ = type(obj)
        return self.fset.__get__(obj, type_)(value)

    def setter(self, func):
        if not isinstance(func, (classmethod, staticmethod)):
            func = classmethod(func)
        self.fset = func
        return self


def classproperty(func):
    if not isinstance(func, (classmethod, staticmethod)):
        func = classmethod(func)

    return ClassPropertyDescriptor(func)

def start_rik():
    """
    Parse command line arguments for RIKPP.

    Returns
    -------
    opt : dict
        A dictionary containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@wkd', 
                        type=str, 
                        default="./", 
                        help="Working directory")
    parser.add_argument('@plot', 
                        action='store_true', 
                        help="Plot the results?")
    parser.add_argument('@tag', 
                        type=str,
                        default='slip', 
                        help="Tag to identify the results")
    parser.add_argument('@L', 
                        type=float,
                        default=15., 
                        help="Fault length [km]")
    parser.add_argument('@W', 
                        type=float,
                        default=15., 
                        help="Fault width [km]")
    parser.add_argument('@nL', 
                        type=int,
                        default=294, 
                        help="Number of grids along fault length")
    parser.add_argument('@nW', 
                        type=int,
                        default=213, 
                        help="Number of grids along fault width")
    parser.add_argument('@nt', 
                        type=int,
                        default=480, 
                        help="Number of time steps")
    parser.add_argument('@dt', 
                        type=float,
                        default=0.025, 
                        help="Time step [s]")
    parser.add_argument('@M0',
                        type=float,
                        default=5.0e+20,
                        help="Seismic moment magnitude [Nm]")
    parser.add_argument('@hL', 
                        type=float,
                        default=12.5, 
                        help="Hypocenter along-length coordinate [km]")
    parser.add_argument('@hW',
                        type=float,
                        default=0.0, 
                        help="Hypocenter along-width coordinate [km]")
    parser.add_argument('@hE',
                        type=float,
                        default=None,
                        help="Hypocenter east coordinate [m]")
    parser.add_argument('@hN',
                        type=float,
                        default=None,
                        help="Hypocenter north coordinate [m]")
    parser.add_argument('@hLON', 
                        type=float, 
                        default=0.0,
                        help="Hypocenter longitude [deg]")
    parser.add_argument('@hLAT', 
                        type=float, 
                        default=0.0,
                        help="Hypocenter latitude [deg]")
    parser.add_argument('@hZ',
                        type=float,
                        default=0.0,
                        help="Hypocenter elevation coordinate [m]")
    parser.add_argument('@strike',
                        type=float,
                        nargs='*',
                        default=[45.0],
                        help="Strike angle(s)")
    parser.add_argument('@dip',
                        type=float,
                        nargs='*',
                        default=[60.0],
                        help="Dip angle(s)")
    parser.add_argument('@rake',
                        type=float,
                        nargs='*',
                        default=[108.0],
                        help="Rake angle(s)")
    parser.add_argument('@sf',
                        type=str,
                        default="slipdistribution.dat",
                        help="Filename of slip distributions")
    parser.add_argument('@fg',
                        type=str,
                        default="xxx.png",
                        help="Figure file name")
    parser.add_argument('@mf',
                        type=str,
                        default='napa2014_moment.h5',
                        help='Moment rate tensor file')
    parser.add_argument('@kf',
                        type=str,
                        default='napa2014_kine.h5',
                        help='Fault parameters file')
    parser.add_argument('@ct',
                        nargs='*',
                        default = ['moment'],
                        help='Fault parameters file')
    parser.add_argument('@snapshots',
                        type=float,
                        default=10.0,
                        help='Time interval between snapshots')
    parser.add_argument('@segL', 
                        type=float, 
                        nargs='*', 
                        default=[0.0],
                        help='Fault segments left-lower corner along-width coordinate')
    parser.add_argument('@segW', 
                        type=float, 
                        nargs='*', 
                        default=[0.0],
                        help='Fault segments left-lower corner along-dip (bottom-up) coordinate')
    parser.set_defaults(plot=True)
    opt = parser.parse_args().__dict__
    
    if (not opt['hE'] or not opt['hN']):
        opt["hE"], opt["hN"], _, _ = utm.from_latlon(opt["hLAT"], opt["hLON"])

    return opt

def get_rotation_tensor(strike: float, dip: float) -> np.float64:
    """
    Compute Rotation tensor from [N,E,D] et to [E,N,Z]
    """
    MatMesh = np.zeros((3,3))

    # See e.g. Aki & Richards, 2002, eq. 3.14
    MatMesh[0,0] =  np.sin(strike)
    MatMesh[0,1] = -np.cos(strike)*np.cos(dip)
    MatMesh[1,0] =  np.cos(strike)
    MatMesh[1,1] =  np.sin(strike)*np.cos(dip)
    MatMesh[1,2] = -np.sin(strike)*np.sin(dip)
    MatMesh[2,1] =  np.sin(dip)
    MatMesh[2,2] =  np.cos(dip)

    # Remove small values
    MatMesh[np.where(MatMesh<=1e-5)]=0.0
    # for i in np.arange(3):
    #     for j in np.arange(3):
    #         if abs(MatMesh[i,j]) < 1e-6:
    #             MatMesh[i,j] = 0.0
    return MatMesh

def moment_computation(M0, time, f, ts, gamma):
    """
    Compute the moment of an earthquake as a function of time.

    Parameters
    ----------
    M0 : float
        The seismic moment at t=0.
    time : np.ndarray
        Array of time values.
    f : float
        The corner frequency of the moment time function.
    ts : float
        The start time of the moment time function.
    gamma : float
        The exponent of the moment time function.

    Returns
    -------
    np.ndarray
        The moment at each time step.
    """
    T 	   = 1.0/f
    moment = np.zeros(len(time))
    for i in np.arange(len(time)):
        t = time[i]
        # Check if the time is greater or equal to the start time
        if (t >= ts):
            # Compute the dimensionless time s
            s = ((t-ts)/T)**gamma
            # Compute the moment at this time step
            moment[i] = M0* (1.0- (1.0+s)*math.e**(-s))
    return moment


def compute_seismic_moment_vectors(strike: float, 
                                   dip: float, 
                                   rake: float) -> tuple[np.float64]:
    """
    Compute the unit-norm normal vector on the fault plane, 
    the unit-norm slip vector in the ENU-xyz reference frame, 
    and the moment matrix.

    Parameters
    ----------
    strike : float
        Strike angle in radians.
    dip : float
        Dip angle in radians.
    rake : float
        Rake angle in radians.

    Returns
    -------
    nv : np.ndarray
        Unit-norm normal vector on the fault plane.
    dv : np.ndarray
        Unit-norm slip vector in the ENU-xyz reference frame.
    M : np.ndarray
        Moment matrix.
    """
    # Compute the unit-norm normal vector on the fault plane
    # x @--> EST
    # y @--> NORD
    # z @--> UP
    nv = np.array([+np.sin(dip)*np.cos(strike),\
                   -np.sin(dip)*np.sin(strike),\
                   +np.cos(dip)])
    nv = nv/np.sqrt(np.dot(nv,nv)) # normalize

    # Compute the unit-norm slip vector in the ENU-xyz reference frame
    dv = np.array([-np.sin(rake)*np.cos(dip)*np.cos(strike) +\
                    np.cos(rake)*np.sin(strike),\
                   +np.sin(rake)*np.cos(dip)*np.sin(strike) +\
                    np.cos(rake)*np.cos(strike),\
                   +np.sin(rake)*np.sin(dip)])
    dv = dv/np.sqrt(np.dot(dv,dv)) # normalize

    # Compute the moment matrix
    M = np.tensordot(nv, dv, axes=0) + np.tensordot(dv, nv, axes=0)

    print('Moment matrix: ')
    print(M[0,:])
    print(M[1,:])
    print(M[2,:])
    return nv, dv, M


def ConvertLine(x): return float(x.rsplit()[1])
def parseRIKrates(filename: str, dimensions: tuple) -> np.float64:
    """Read moment rate file from RIKsrf

    Arguments:
        filename (string): rate filename from RIKsrf
        dimensions (tuple): nL, nW, nt
    Returns:
        rate (numpy.array): rate array 
    """
    nL, nW, nt = dimensions
    rate = np.zeros((nL*nW, nt), 
                    dtype=np.float64)
    fid = open(filename, 'r')
    FileContent = fid.readlines()
    
    rate = np.vstack([np.array(list(map(ConvertLine,
                                        FileContent[i*(nt+2):
                                            (i+1)*(nt+2)-2]))) 
                    for i in range(nL*nW)]
                    ).reshape(dimensions, order='F')
    return rate

def make_colormap(seq):
    """Return a LinearSegmentedColormap

    seq is a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).

    The output colormap is a LinearSegmentedColormap with the colors specified
    in the input sequence.

    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            # item is a float, so we need to interpolate the color
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
        else:
            # item is an RGB tuple, so we add it directly to the cdict
            cdict['red'].append([item[0], item[0], item[0]])
            cdict['green'].append([item[1], item[1], item[1]])
            cdict['blue'].append([item[2], item[2], item[2]])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def readhypo(filename):
    """
    Reads a hypocenter file and returns the coordinates of the hypocenter.

    Parameters
    ----------
    filename : str
        The name of the file containing the hypocenter coordinates

    Returns
    -------
    x, y : numpy arrays
        The x and y coordinates of the hypocenter

    """
    x = np.genfromtxt(filename, usecols=0)
    y = np.genfromtxt(filename, usecols=1)
    return x, y

def read_subsources (filename):
    """
    Reads a file containing the coordinates of the subsources.

    Parameters
    ----------
    filename : str
        The name of the file containing the coordinates of the subsources

    Returns
    -------
    xgrid, ygrid : numpy arrays
        The x and y coordinates of the subsources
    ssources : numpy array
        The subsources coordinates

    """
    f = open (filename, 'r')
    xgrid     = np.genfromtxt(filename, usecols=0)
    ygrid     = np.genfromtxt(filename, usecols=1)
    ssources  = np.genfromtxt(filename, usecols=2)
    return xgrid, ygrid, ssources


def readfileslip (NSR, filename):
    """
    Reads a file containing the slip distribution over the fault plane.

    Parameters
    ----------
    NSR : int
        The number of points where the slip is evaluated
    filename : str
        The name of the file containing the slip distribution

    Returns
    -------
    xgrid, ygrid : numpy arrays
        The x and y coordinates of the points where the slip is evaluated
    slip : numpy array
        The slip distribution

    """
    xgrid = np.zeros((NSR))
    ygrid = np.zeros((NSR))
    slip  = np.zeros((NSR))
    f = open (filename, 'r')
    for i in np.arange(NSR):
        string   = f.readline()
        degerler = [float(j) for j in string.split()]
        xgrid[i] = degerler[0]
        ygrid[i] = degerler[1]
        slip [i] = degerler[2]
    return xgrid, ygrid, slip

def plotslip(nrows,ncols,LF,WF,slip,kaydet,figname,nucx,nucy):
    """
    Plot the maximum slip distribution over the fault plane.

    Parameters
    ----------
    nrows, ncols : int
        The number of points where the slip is evaluated
    LF, WF : float
        The length and width of the fault
    slip : numpy array
        The maximum slip distribution
    kaydet : bool
        If True, the plot will be saved as a file
    figname : str
        The name of the file
    nucx, nucy : float
        The x and y coordinates of the nucleation point

    """
    print('\n\n')
    print('Plotting max-slip distribution')
    print(nrows,ncols,LF,WF)
    # Making a colormap
#    c    = mcolors.ColorConverter().to_rgb
#    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
#          '#fcab6c', '#f50000']
#    cc1 = np.linspace(0,1,6)
#    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0],
#                          c(cc[1]), c(cc[2]), cc1[1],
#                          c(cc[2]), c(cc[3]), cc1[2],
#                          c(cc[3]), c(cc[4]), cc1[3],
#                          c(cc[4]), c(cc[5]), cc1[4],
#                          c(cc[5]), c(cc[6]), cc1[5],
#                          c(cc[6])])
    sns.set_style('whitegrid')
    cmap = cm.RdBu_r

    fig = p.figure(figsize=(18,10))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$x$ Along strike [km]', fontsize=17)
    ax.set_ylabel(r'$y$ Along up-dip [km]', fontsize=17)
    
    ax.set_xlim([0.,LF])
    ax.set_ylim([0.,WF])
    vmin = slip.min()
    vmax = slip.max()
    print('Slip: min {:>.1f} - max {:>.1f}\n'.format(vmin,vmax))
    grid = slip.reshape((nrows, ncols))
    im = plt.imshow(grid, extent=(0.0,LF,0.0,WF),interpolation='bilinear',cmap=cmap,origin='lower')
    formatter = FormatStrFormatter('%.1f') #LogFormatter(10,labelOnlyBase=False)
    #norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cb = fig.colorbar(im,shrink=0.5,aspect=10,pad=0.01,
            ticks=np.linspace(vmin,vmax,11),format=formatter)
    cb.set_label(r'$\Delta u$ [m]',labelpad=20,y=0.5,rotation=90,fontsize=17)
    p.setp(cb.ax.yaxis.get_ticklabels(), fontsize=16)
    ax.plot(nucx,nucy, marker='*',color='red',markersize=20)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    topname  = r"$\max\Delta u$ = {:>.2f} m".format(max(slip))
    plt.title(topname,fontsize=20)
    if kaydet:
        print("{:>s}_slipmax.png".format(figname))
        plt.savefig("{:>s}_slipmax.png".format(figname),dpi=500,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def plotsubsources(x, y, radius, LF, WF, kaydet, figname, nucx, nucy):
    """
    Plot the sub-sources on the fault plane
    
    Parameters
    ----------
    x, y : numpy arrays
        The coordinates of the sub-sources
    radius : numpy array
        The radius of each sub-source
    LF, WF : float
        The length and width of the fault plane
    kaydet : bool
        If True, the plot will be saved as a file
    figname : str
        The name of the file
    nucx, nucy : float
        The x and y coordinates of the nucleation point
    
    """
    print('Plotting sub-sources\n')
    
    # Making a colormap
    c    = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
          '#fcab6c', '#f50000']
    
    #cc1 = np.logspace(np.log10(0.25),np.log10(0.95),6)
    cc1 = np.linspace(0,1,6)
    
    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0], 
                          c(cc[1]), c(cc[2]), cc1[1],
                          c(cc[2]), c(cc[3]), cc1[2],
                          c(cc[3]), c(cc[4]), cc1[3],
                          c(cc[4]), c(cc[5]), cc1[4],
                          c(cc[5]), c(cc[6]), cc1[5],
                          c(cc[6])])
    
    # Figure parameters
    sns.set_style('whitegrid')
    fig = p.figure(figsize=(18,10))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Along strike [km]', fontsize=17)
    ax.set_ylabel('Along up-dip [km]', fontsize=17)
    ax.set_xlim([0,1.1*LF])
    ax.set_ylim([-0.2,1.1*WF])
    
    # Plot the sub-sources
    i = 0
    for xx,yy in zip(x,y):
        # Create a circle for each sub-source
        circ = Circle((xx,yy),radius[i], fill=False, lw=1)
        ax.add_patch(circ)	   
        i = i+1
    
    # Number of sub-sources
    Nsub = len(radius)
    topname  = 'Total number of subsources = '+ '%d' % Nsub
    plt.title(topname+'\n', fontsize=20)
    ax.plot(nucx,nucy, marker='*', color='red',markersize=20)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    
    # Save the plot if needed
    if kaydet:
        fig.savefig(figname, dpi=300)
    else:
        plt.show()

def plotgrid(nrows,ncols,LF,WF,kaydet,figname,hypox,hypoy,x,y):
    """
    Plot the 2D grid of the fault plane.

    Parameters
    ----------
    nrows, ncols : int
        The number of points where the slip is evaluated
    LF, WF : float
        The length and width of the fault
    kaydet : bool
        If True, the plot will be saved as a file
    figname : str
        The name of the file
    hypox, hypoy : float
        The x and y coordinates of the nucleation point
    x, y : numpy arrays
        The coordinates of the sub-sources

    """
    print('Creating the 2D grid...')

    # Create a colormap
    c    = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
          '#fcab6c', '#f50000']
    #cc1 = np.logspace(np.log10(0.25),np.log10(0.95),6)
    cc1 = np.linspace(0,1,6)
    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0], 
                          c(cc[1]), c(cc[2]), cc1[1],
                          c(cc[2]), c(cc[3]), cc1[2],
                          c(cc[3]), c(cc[4]), cc1[3],
                          c(cc[4]), c(cc[5]), cc1[4],
                          c(cc[5]), c(cc[6]), cc1[5],
                          c(cc[6])])
    # Figure parameters
    sns.set_style('whitegrid')
    fig = p.figure(figsize=(18,10))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Along strike [km]', fontsize=17)
    ax.set_ylabel('Along up-dip [km]', fontsize=17)
    ax.set_xlim([0,LF])
    ax.set_ylim([-0.2,WF])
    plt.scatter(x,y,s=2,c='gray')
    ax.plot(hypox,hypoy, marker='*', color='red', markersize=25)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    topname  = '2D grid for RIK model'
    plt.title(topname+'\n', fontsize=20)
    if kaydet:
        fig.savefig(figname, dpi=300)
    else:
        plt.show()
