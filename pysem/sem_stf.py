# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Tool to produce different wavelets
"""
# Requested modules
import matplotlib as mpl
#mpl.use('Agg') # de-comment if pyplot raise errors
from matplotlib import pyplot as plt
import argparse
import numpy as np
from scipy.signal import ricker
import warnings


# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

"""
Find 2^n that is equal to or greater than.
"""
def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n

def plot_write(stf,sff,tag):
    np.savetxt(tag+'_stf.txt',stf,fmt='%.15f', delimiter=' ')
    np.savetxt(tag+'_sff.txt',sff.view(float)[:,[0]+list(range(2,4))],fmt='%.15f', delimiter=' ')
    plt.figure(figsize=(12,5))
    plt.plot(stf[:,0],stf[:,1],linewidth=3,label=tag)
    plt.xlabel(r'$\mathbf{t [s]}$',fontsize=14)
    plt.ylabel(r'$\mathbf{STF}$',fontsize=14)
    plt.legend()
    plt.savefig(tag+'_stf.png',bbox_inches='tight')
    plt.close()
    plt.figure(figsize=(5,5))
    plt.loglog(np.real(sff[:,0]),np.abs(sff[:,1]),linewidth=3,label=tag)
    plt.xlabel(r'$\mathbf{f [Hz]}$',fontsize=14)
    plt.ylabel(r'$\mathbf{STF}$',fontsize=14)
    plt.xlim(10**-1,10**2)
    plt.legend()
    plt.savefig(tag+'_sff.png',bbox_inches='tight')
    plt.close()
    
def ricker(vtm,ts,wd,tag=None,plot=False):
    dtm = vtm[1]-vtm[0]
    npt = vtm.size
    nfr = nextpow2(npt)
    vfr = np.linspace(0,1./dtm,nfr).reshape(-1,1)
    
    sig = (vtm-ts*np.ones_like(vtm))**2/wd
    stf = 2.*(1.+2.*sig)*np.exp(sig)/wd
    stf[vtm<=ts] = 0.
    
    sff = np.fft.fft(stf,n=nfr,axis=0)
    stf = np.concatenate((vtm,stf),axis=-1)
    sff = np.concatenate((vfr,sff),axis=-1)
    
    if plot:
        plot_write(stf,sff,tag)
        
    return stf,sff
    
def gaussian(vtm,ts,wd,tag=None,plot=False):
    dtm = vtm[1]-vtm[0]
    npt = vtm.size
    nfr = nextpow2(npt)
    vfr = np.linspace(0,1./dtm,nfr).reshape(-1,1)
    
    sig = (vtm-ts*np.ones_like(vtm))/wd
    
    stf = np.exp(-sig**2)
    
    sff = np.fft.fft(stf,n=nfr,axis=0)
    stf = np.concatenate((vtm,stf),axis=-1)
    sff = np.concatenate((vfr,sff),axis=-1)
    
    if plot:
        plot_write(stf,sff,tag)
    return stf,sff

def spice_bench(vtm,ts,wd,k,tag=None,plot=False):
    dtm = vtm[1]-vtm[0]
    npt = vtm.size
    nfr = nextpow2(npt)
    vfr = np.linspace(0,1./dtm,nfr).reshape(-1,1)
       
    sig = (vtm-ts*np.ones_like(vtm))/wd
    print(1./wd,ts)
    
    stf = np.ones_like(vtm)-(np.ones_like(vtm)+sig)*np.exp(-sig)
    stf[vtm<=ts] = 0.
    
    sff = np.fft.fft(stf,n=nfr,axis=0)
    stf = np.concatenate((vtm,stf),axis=-1)
    sff = np.concatenate((vfr,sff),axis=-1)
    
    if plot:
        plot_write(stf,sff,tag)
    
    return stf,sff


if __name__== '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--btm',type=float,default=0.0,help="Begin time step")
    parser.add_argument('--ftm',type=float,default=0.2,help="Final time step")
    parser.add_argument('--npt',type=int,default=10000,help="Number of time steps")
    parser.add_argument('--model',type=str,nargs="*",default=["spice_bench"],help="Wavelet model")
    parser.add_argument('--ts',type=float,nargs="*",default=[0.2],help="Time shift")
    parser.add_argument('--wd',type=float,nargs="*",default=[0.0002],help="Pulse width")
    parser.add_argument('--fc',type=float,nargs="*",default=[35.250],help="Dominant frequency")
    parser.add_argument('--k',type=int,nargs="*",default=[1],help="Dominant frequency")
    parser.add_argument('--ps',action='store_true',default=False,help='Plot single?')
    parser.add_argument('--pc',action='store_true',default=False,help='Plot comparison?')

    opt = parser.parse_args().__dict__
     
    vtm = np.linspace(opt['btm'],opt['ftm'],opt['npt']+1).reshape(-1,1)
    
    stf={}
    sff={}
    for i in range(len(opt['model'])):
        md = opt['model'][i].lower()
        ts = opt['ts'][i]
        fc = opt['fc'][i]
        try:
            k = opt['k'][i]
        except:
            warnings.warn("k not defined")
            pass
        tag = "md={}-ts={:>5.2f}s-fc={:>5.2f}Hz".format(md,round(ts,2),round(fc,2))
        if 'gaussian' in md:
            wd = 1./fc
            stf[tag],sff[tag]=gaussian(vtm,ts,wd,tag,opt['ps'])
        elif 'ricker' in md:
            wd = -1./(np.pi**2*fc**2)
            stf[tag],sff[tag]=ricker(vtm,ts,wd,tag,opt['ps'])
        elif 'spice_bench' in md:
            wd = 1./fc
            stf[tag],sff[tag]=spice_bench(vtm,ts,wd,k,tag,opt['ps'])
    
    if opt['pc']:
        plt.figure(figsize=(12,5))
        for i in range(len(opt['model'])):
            md = opt['model'][i].lower()
            ts = opt['ts'][i]
            fc = opt['fc'][i]
            tag = "md={}-ts={:>5.2f}s-fc={:>5.2f}Hz".format(md,round(ts,2),round(fc,2))
            plt.plot(stf[tag][:,0],stf[tag][:,1],linewidth=3,
                     label=tag)
        plt.xlabel(r'$\mathbf{t [s]}$',fontsize=14)
        plt.ylabel(r'$\mathbf{STF}$',fontsize=14)
        plt.legend()
        plt.savefig('compare_stf.png',bbox_inches='tight')
        plt.close()
        plt.figure(figsize=(5,5))
        for i in range(len(opt['model'])):
            md = opt['model'][i].lower()
            ts = opt['ts'][i]
            fc = opt['fc'][i]
            tag = "md={}-ts={:>5.2f}s-fc={:>5.2f}Hz".format(md,round(ts,2),round(fc,2))
            plt.loglog(np.real(sff[tag][:,0]),np.abs(sff[tag][:,1]),linewidth=3,
                       label=tag)
        plt.xlabel(r'$\mathbf{f [Hz]}$',fontsize=14)
        plt.ylabel(r'$\mathbf{STF}$',fontsize=14)
        plt.xlim(10**-1,10**2)
        plt.legend()
        plt.savefig('compare_sff.png',bbox_inches='tight')
        plt.close()