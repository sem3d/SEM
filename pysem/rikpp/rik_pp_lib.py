# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Set of subroutines to post-process RIK output
"""

#=======================================================================
# Required modules
#=======================================================================
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
from   matplotlib.patches import Polygon
from   matplotlib.collections import PatchCollection
from   matplotlib.colors import LinearSegmentedColormap, Normalize
from   matplotlib.cm import ScalarMappable
from   matplotlib import ticker
from   matplotlib.colors import LogNorm
from   matplotlib.image import NonUniformImage
from   matplotlib.ticker import LogFormatter, FormatStrFormatter
from   matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm,rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter
rc('text', usetex=True)
import seaborn as sns
from   collections import OrderedDict 

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

def start_rik():
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@wkd',type=str,default="./",help="Working directory")
    parser.add_argument('@plot',action='store_true',help="plot?")
    parser.add_argument('@tag',type=str,default='slip',help="tag")
    parser.add_argument('@L',type=float,default=15.,help="fault length [km]")
    parser.add_argument('@W',type=float,default=15.,help="fault width [km]")
    parser.add_argument('@nL',type=int,default=294,help="Number of grids along fault length")
    parser.add_argument('@nW',type=int,default=213,help="Number of grids along fault width")
    parser.add_argument('@nt',type=int,default=480,help="Number of time steps")
    parser.add_argument('@dt',type=float,default=0.025,help="Time step [s]")
    parser.add_argument('@hL',type=float,default=12.5,help="Hypocenter along-length coordinate [km]")
    parser.add_argument('@hW',type=float,default=0.0,help="Hypocenter along-width coordinate [km]")
    parser.add_argument('@hE',type=float,default=0.0,help="Hypocenter east coordinate [m]")
    parser.add_argument('@hN',type=float,default=0.0,help="Hypocenter north coordinate [m]")
    parser.add_argument('@hZ',type=float,default=0.0,help="Hypocenter elevation coordinate [m]")
    parser.add_argument('@strike',type=float,default=45.0,help="Strike")
    parser.add_argument('@dip',type=float,default=60.0,help="Dip")
    parser.add_argument('@rake',type=float,default=108.0,help="Rake")
    parser.add_argument('@sf',type=str,default="slipdistribution.dat",help="Filename of slip distributions")
    parser.add_argument('@fg',type=str,default="xxx.png",help="Figure file name")
    parser.add_argument('@mf',type=str,default='napa2014_moment.h5',help='Moment rate tensor file')
    parser.add_argument('@kf',type=str,default='napa2014_kine.h5',help='Fault parameters file')
    parser.add_argument('@ct',nargs='*',default = ['moment'],help='Fault parameters file')
    parser.add_argument('@snapshots',type=float,default=10.,help='Time interval between snapshots')
    parser.set_defaults(plot=True)
    opt = parser.parse_args().__dict__
    return opt

def get_rotation_tensor(aS,aD):
    r'''From [N,E,D] et to [E,N,Z]
    '''
    MatMesh = np.zeros((3,3))
    MatMesh[0,0] = -np.cos(aS)* np.cos(aD)
    MatMesh[0,1] = +np.sin(aS)
    MatMesh[1,0] = +np.sin(aS)* np.cos(aD)
    MatMesh[1,1] = +np.cos(aS)
    MatMesh[2,2] =  -1
    for i in np.arange(3):
        for j in np.arange(3):
            if abs(MatMesh[i,j]) < 1e-15:
                MatMesh[i,j] = 0.0
    return MatMesh

def moment_computation(M0,time,f,ts,gamma):
    T 	   = 1.0/f
    moment = np.zeros(len(time))
    for i in np.arange(len(time)):
        t = time[i]
        if (t >= ts):
            s = ((t-ts)/T)**gamma
            moment[i] = M0* (1.0- (1.0+s)*math.e**(-s))
    return moment

def yeni_data_seti (data,dim,chunk, attrib,output):
    ''' Attribuer un set de data de dimension dim
    comme attribut de attrib dans un fichier output de hdf5 '''
    dataset = output.create_dataset(attrib, dim, chunks=chunk)
    dataset = data

def vecteurs(aD,aS,aR,output):
    ''' Calcul du vecteur normal Vnormal sur le plan de faille, et du vecteur unitaire
    de glissement Vslip dans le repere de reference et les attribuer dans le fichier 
    output de hdf5 '''
    # En supposant que (puisque l'on multiplie avec MatMesh)
    # x @--> EST
    # y @--> NORD
    # z @--> UP
    Vnormal = np.array([+np.sin(aD)*np.cos(aS),\
                        -np.sin(aD)*np.sin(aS),\
                        +np.cos(aD)])
    Vslip   = np.array([-np.sin(aR)*np.cos(aD)*np.cos(aS) + np.cos(aR)*np.sin(aS),\
                        +np.sin(aR)*np.cos(aD)*np.sin(aS) + np.cos(aR)*np.cos(aS),\
                        +np.sin(aR)*np.sin(aD)])
    output.attrs['Vnormal'] = Vnormal
    output.attrs['Vslip']   = Vslip
    print("Vector normal to the fault : {}".format(Vnormal))
    print("Vector of the slip         : {}".format(Vslip))
    M = np.zeros((3,3))
    for i in np.arange(3):
        for j in np.arange(3):
            M[i,j] = Vnormal[i]*Vslip[j]+Vnormal[j]*Vslip[i]
    print('Moment matrix of SEM3D code: ')
    print(M[0,:])
    print(M[1,:])
    print(M[2,:])

def read_moment_rate_RIK(dosya,dim):
    ''' Lire la vitesse de moment par un fichier dosya et 
    l'assigner a l'array de dimension dim '''
    
    Mrate = np.zeros((dim[0]*dim[1], dim[2]))
    f = open (dosya, 'r')
    # Pour tous les points
    for i in np.arange(dim[0] * dim[1]):
        # Histoire temporelle
        for t in np.arange(dim[2]):
            string = f.readline()
            dagit = string.rsplit()
            Mrate[i, t] = float(dagit[1])
        f.readline(); f.readline()
    Mratenew = np.reshape(Mrate, dim, order='F')
    return Mratenew

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def readhypo(filename):
    x = np.genfromtxt(filename,usecols=0)
    y = np.genfromtxt(filename,usecols=1)
    return x,y

def read_subsources (filename):

    f = open (filename, 'r')
    xgrid     = np.genfromtxt(filename, usecols=0)
    ygrid     = np.genfromtxt(filename, usecols=1)
    ssources  = np.genfromtxt(filename, usecols=2)
    return xgrid, ygrid, ssources


def readfileslip (NSR, filename):
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
    topname  = r"$\max\Delta u$ [m] = {{{:>.2f}}}".format(max(slip))
    plt.title(topname,fontsize=20)
    if kaydet:
        print("{:>s}_slipmax.png".format(figname))
        plt.savefig("{:>s}_slipmax.png".format(figname),dpi=500,bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def plotsubsources(x,y,radius,LF,WF,kaydet,figname,nucx,nucy):

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
    i = 0
    for xx,yy in zip(x,y):
        circ = Circle((xx,yy),radius[i], fill=False, lw=1)
        ax.add_patch(circ)	   
        i = i+1
    Nsub = len(radius)
    topname  = 'Total number of subsources = '+ '%d' % Nsub
    plt.title(topname+'\n', fontsize=20)
    ax.plot(nucx,nucy, marker='*', color='red',markersize=20)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    if kaydet:
        fig.savefig(figname, dpi=300)
    else:
        plt.show()

def plotgrid(nrows,ncols,LF,WF,kaydet,figname,hypox,hypoy,x,y):

    print('Creating the 2D grid...')
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
