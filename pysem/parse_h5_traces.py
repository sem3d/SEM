# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to parse and post-process SEM3D traces
    
    Ex. : Parse hdf5 traces called Uobs, for all directions x y z, monitor 0 and 1 (in stations.txt) 
        and plot them
        
        python3 ParseSEM3DH5Traces_h5_traces.py @@wkd ./traces/ @@fmt h5 @@nam Uobs @@var Displ @@rdr x y z @@mon 0 1 2 @@plt
"""

# Required modules
import matplotlib as mpl
mpl.use('Agg') # de-comment if pyplot raise errors
from matplotlib import pyplot as plt
import argparse
import os.path as osp
from os.path import join as osj
from collections import OrderedDict
import numpy as np
import h5py
from scipy.interpolate import interp1d
from scipy.signal import kaiserord, firwin, lfilter, convolve, medfilt, decimate
import glob
import copy

cmp_dict = [{'p':0}]
cmp_dict.append({'x':0,'y':1,'z':2})
cmp_dict.append({'xx':0,'yy':1,'zz':2,'xy':3,'xz':4,'yz':5})
components = {
    'Pressure':cmp_dict[0],
    'Displ':cmp_dict[1],
    'Veloc':cmp_dict[1],
    'Accel':cmp_dict[1],
    'StressDev':cmp_dict[2],
    'EpsDev':cmp_dict[2]
    }
units = {'Displ':'m','Veloc':'m/s','Accel':'m/s^2',
         'Pressure':'Pa','StressDev':'Pa','EpsDev':'1','EpsDevPl':'1'}

class SEM3DMonitor(object):
    def __init__(self, name, fmt, data={}, \
                 nt = 0,nc = 0, dTime=0.,\
                 var=OrderedDict({}), \
                 var_avl =OrderedDict({}), \
                 comp = ['x','y','z']):
        self.name = name
        self.fmt  = fmt
        self.data = data
        self.var  = var
        self.var_avl  = var_avl
        self.comp = comp
        self.pos = {'0':np.array([0.,0.,0.])}
        self.nc = nc
        self.nt = nt
        self.dTime = dTime
        
    def __call__(self, **kwargs):    
        self.__dict__.update(kwargs)
    
    def c2d(self,k,c):
        v = self.var_avl[k]
        c2d = [v[i] for i in self.comp[k]]
        return c2d.index(c)
    
    def add_capteur(self,data,c):
        for k,v in self.var.items():
            if k != 'Time':
                self.data[k][:,:,c] = data[:,[v[i]+1 for i in self.comp[k]]]
                
    def add_position(self,data,c):
        self.set_coord(c,data[0],data[1],data[2])
            
    def set_coord(self,c,x,y,z):
        self.pos[str(c)] = np.array([x,y,z])
        
    def set_dimensions(self):
        self.data = {}
        self.position = {}
        for k in self.var.keys():
            if k!='Time':
                self.data[k] = np.empty((self.nt,len(self.comp[k]),self.nc),dtype=np.float_)
        for c in range(self.nc):
            self.position[str(c)] = np.empty((3,),dtype=np.float_)         
                
    def set_time(self,Time):
        self.Time = Time
        self.set_dTime()
        self.set_nt()
        
    def set_nt(self):
        if self.fmt=='h5' and 'Time' in self.var:
            try:
                self.nt = self.Time.size
            except:
                raise('Time vector not ParseSEM3DH5Tracesd!')
            
    def set_dTime(self):
        if self.fmt=='h5' and 'Time' in self.var.keys():
            try:
                self.dTime = self.Time[-1]-self.Time[-2]
            except:
                raise('Time vector not ParseSEM3DH5Tracesd!')
            
        else:
            # on ne prend pas [1]-[0] car certains filtres decalent le T0
            N=self.data.shape[0]/2
            return self.data[N+1,0] - self.data[N,0]
    
    def total_stress(self,rdr,cpt='all'):
        if cpt=='all':
            cpt = [c for c in range(self.nc)]
        tts = OrderedDict({})
        for m in cpt:
            ts = {}
            for cc in rdr:
                c = components['StressDev'][cc]
                if c in self.comp['StressDev']:
                    if c<3:
                        ts[cc] = self.data['Pressure'][:,0,m]+\
                            self.data['StressDev'][:,self.c2d('StressDev',c),m]
                    else:
                        ts[cc] = self.data['StressDev'][:,self.c2d('StressDev',c),m]
            tts[str(m)] = ts
        return tts

    def plot(self,var,rdr,mon,hfg=None,svf=False,**kwargs):
        for m in mon:
            for v in var:
                if v in self.var.keys():
                    for c in rdr:
                        c = components[v][c]
                        if c in self.comp[v]: 
                            if hfg:
                                plt.figure(hfg.number)
                            else:
                                plt.figure(figsize=[10,5])
                            plt.plot(self.Time,self.data[v][:,self.c2d(v,c),m])
                            plt.xlim(self.Time[0],self.Time[-1])
                            plt.xlabel(r'$\mathbf{t [s]}$',fontsize=14)
                            plt.ylabel(r'$\mathbf{{ {{{vv}}} [{{{uu}}}] }}$'.format(vv=v,uu=units[v]),fontsize=14)
                            if svf:
                                plt.savefig(self.name+'_'+str(m)+'_'+str(v)+'_'+str(c)+'.png',\
                                            dpi=300,bbox_inches='tight')
                                plt.close()
                            else:
                                return plt.gcf()
                        else:
                            raise ValueError('Component '+c+' not ParseSEM3DH5Tracesd!')
                else:
                    raise ValueError('Variable '+v+' not ParseSEM3DH5Tracesd!')
            

    def filt_low_pass(self, fmax, fband=2.0):
        sample_rate = 1./self.dt()

        # The Nyquist rate of the signal.
        nyq_rate = sample_rate / 2.0

        # The desired width of the transition from pass to stop,
        # relative to the Nyquist rate.  We'll design the filter
        # with a 5 Hz transition width.
        width = fband/nyq_rate

        # The desired attenuation in the stop band, in dB.
        ripple_db = 60.0

        # Compute the order and Kaiser parameter for the FIR filter.
        N, beta = kaiserord(ripple_db, width)

        # The cutoff frequency of the filter.
        cutoff_hz = fmax

        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        taps = firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))

        # Use lfilter to filter x with the FIR filter.
        data2 = np.array(self.data, copy=True)
        for comp in 1,2,3:
            data2[:,comp] = lfilter(taps, 1.0, self.data[:,comp])

        return self.filt_capteur(data2, "LP(%f Hz)" % fmax)

    def filt_capteur(self, data, pfx):
        cpt = self.copy()
        return cpt

    def filt_ma(self, n, decimate=True):
        nn = 2*n+1
        flt = np.ones((nn,))/float(nn)
        data = np.array(self.data, copy=True)
        for comp in 1,2,3:
            data[:,comp] = convolve(self.data[:,comp], flt, mode="same")
        if decimate:
            data = np.array(data[n::nn,:], copy=True)
        return self.filt_capteur(data, "MA(%d)" % nn)

    def filt_median(self, n, decimate=True):
        data = np.array(self.data, copy=True)
        nn = 2*n+1
        for comp in 1,2,3:
            data[:,comp] = medfilt(self.data[:,comp], nn)
        if decimate:
            data = np.array(data[n::nn,:], copy=True)
        return self.filt_capteur(data, "MED(%d)" % nn)

    def filt_decimate(self, fac, scl):
        cpt = copy.deepcopy(self)
        cpt.name = cpt.name+'_dec'
        
        dtm_rsmpl = cpt.dTime/fac
        ntm_rsmpl = int(cpt.nt*fac)
        vtm_rsmpl = dtm_rsmpl*np.arange(0,ntm_rsmpl)
        cpt.set_time(vtm_rsmpl)
        for v in self.var.keys():
            if v != 'Time':
                f = interp1d(self.Time,self.data[v],axis=0)
                cpt.data[v] = f(cpt.Time)
        return cpt

def ParseSEM3DH5Traces(wkd='./',fmt='h5',var=[''],rdr=['x','y','z'],nam='all',**kwargs):    
    if fmt == 'h5':
        fls = glob.glob(osj(wkd,'*.h5'))
        fn = fls[0]
        f = h5py.File(fn,"r+")
        tmp = [tuple(a.split()) for a in f['Variables'][...].tolist()]
        var_avl = OrderedDict({})
        var_ok = OrderedDict({'Time':0})
        for vi in range(len(tmp)):
            v = tmp[vi]
            v = tuple([vb.decode("utf-8") for vb in v])
            if len(v)>2:
                v=(v[0]+v[1],v[-1])
            if v[0] not in var_avl.keys():
                var_avl[v[0]] = [int(v[-1])-1]
            else:
                var_avl[v[0]].append(int(v[-1])-1)
            if v[0] in var:
                if v[0] not in var_ok.keys():
                    var_ok[v[0]] = [vi-1]
                else:
                    var_ok[v[0]].append(vi-1)

        # define components for each variable
        cmp_ok = {}
        for k in var_ok.keys():
            if k!='Time':
                cmp_ok[k] = [components[k][c] for c in rdr if c in components[k]]
                cmp_ok[k].sort()
        
        if 'all' in nam:
            nam = [c for c in nam]
        cpts = {}
        for n in nam:
            cpt = SEM3DMonitor(name=n,fmt=fmt,var=var_ok,var_avl=var_avl,comp=cmp_ok)
            flag=False
            for ds in f.items():
                if 'Variables' not in ds[0]:
                    if 'pos' not in ds[0] :
                        # ParseSEM3DH5Traces time vector and construct time properties
                        if not flag: 
                            cpt.set_time(ds[1][...][:,0])
                            flag=True
            f.close()
            # ParseSEM3DH5Traces h5 files to find total number of capteurs
            idx = {}
            nc = 0
            for fn in fls:
                idx[fn] = []
                f = h5py.File(fn,"r+")
                for ds in f.items():
                    if 'Variables' not in ds[0] and (n in ds[0] or 'all' in n):
                        if 'pos' in ds[0]:
                            nc += 1
                            idx[fn].append(int(ds[0].split('_')[-2]))
                f.close()
            cpt(nc=nc,idx=idx)
            cpt.set_dimensions()
        
            # loop over h5 files
            for fn in fls:
                idx[fn] = []
                f = h5py.File(fn,"r+")
                for ds in f.items():
                    if 'Variables' not in ds[0] and (n in ds[0] or 'all' in n):
                        if 'pos' in ds[0]:
                            cpt.add_position(ds[1],int(ds[0].split('_')[-2]))
                        else:
                            cpt.add_capteur(ds[1],int(ds[0].split('_')[-1]))
                f.close()
            cpts[n]=cpt
        return cpts

def ExportPickle(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, -1)# pickle.HIGHEST_PROTOCOL)
    output.close()

def ParseOptions():
    OptionParser = argparse.ArgumentParser(prefix_chars='@')
    OptionParser.add_argument('@w','@@wkd',type=str,default='../test/traces/',help='Database main directory')
    OptionParser.add_argument('@f','@@fmt',type=str,default='h5',help='Database format')
    OptionParser.add_argument('@n','@@nam',type=str,nargs='+',default=['all'],help = 'Name of the set(s) of monitors')
    OptionParser.add_argument('@v','@@var',type=str,nargs='+',default=['Displ','Veloc','Accel'],help='Output variables')
    OptionParser.add_argument('@r','@@rdr',type=str,nargs='+',default=['x','y','z'],help='Motion components')
    OptionParser.add_argument('@m','@@mon',type=int,nargs='+',default=[0],help='SEM3DMonitor number')
    OptionParser.add_argument('@p','@@plt',action='store_true',default=True,help='Plot?')
    options = OptionParser.parse_args().__dict__
    return options

if __name__=='__main__':
    
    options = ParseOptions()
    
    print("Parse {} database ({}*.{})\nVariables: {}-Comp: {}".format(options['nam'],options['wkd'],options['fmt'],
                                                              options['var'],options['rdr']))
    
    # ParseSEM3DH5Traces database
    stream = ParseSEM3DH5Traces(**options)
    print("Database ParseSEM3DH5Traces!\n")

    # Plot traces
    if options['plt']:
        if 'all' in stream:
            stream['all'].plot(**options,svf='plot.png')
        else:
            for n in options['nam']:
                stream[n].plot(**options,svf='plot.png')
    
    # Compute total stresses if present
    try:
        ts = sem_db.total_stress(rdr=['xx','yy','zz','xy','xz','yz'])
    except:
        print("Warning: no stress field defined in database")
        pass
    
