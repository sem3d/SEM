# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to parse and post-process SEM3D traces
    
    Ex. : Parse hdf5 traces called Uobs, for all directions x y z, monitor 0 and 1 (in stations.txt) 
        and plot them
        
        python3 parse_sem3d_traces.py @@wkdir ./traces/ @@format h5 @@names Uobs @@variables Displ @@components x y z @@monitors 0 1 2 @@plt
"""

# Required modules
import sys
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
from tqdm import tqdm


# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2023, LMPS UMR CNRS 9026, CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

cmp_dict = [{'p':0}]
cmp_dict.append({'x':0,'y':1,'z':2})
cmp_dict.append({'xx':0,'yy':1,'zz':2,'xy':3,'xz':4,'yz':5})
components_dict = {
    'Pressure':cmp_dict[0],
    'Displ':cmp_dict[1],
    'Veloc':cmp_dict[1],
    'Accel':cmp_dict[1],
    'StressDev':cmp_dict[2],
    'EpsDev':cmp_dict[2],
    }
units = {'Displ':'m','Veloc':'m/s','Accel':'m/s^2',
         'Pressure':'Pa','StressDev':'Pa','EpsDev':'1','EpsDevPl':'1'}

class SEM3DMonitor(object):
    def __init__(self, 
                 name, 
                 format, 
                 data={}, 
                 position = {},
                 nt = 0, nc = 0, dTime=0.,
                 variables=OrderedDict({}), 
                 variables_available =OrderedDict({}),
                 components = ['x','y','z']):
        """
        Initialize SEM3DMonitor object.

        Parameters
        ----------
        name : str
            Name of SEM3DMonitor object.
        format : str
            Format of data. Only 'h5' is supported.
        data : dict
            Dictionary of data arrays.
        nt : int
            Number of time steps.
        nc : int
            Number of capteurs.
        dTime : float
            Time step.
        variables : OrderedDict
            Dictionary of variables.
        variables_available : OrderedDict
            Dictionary of available variables.
        components : list of str
            List of components for each variable.
        """
        self.name = name
        self.format  = format
        self.data = data
        self.pos = position
        self.variables = variables
        self.variables_available = variables_available
        self.components = components
        self.nc = nc
        self.nt = nt
        self.dTime = dTime
        self.TimeSet = False
        
    def __call__(self, **kwargs):    
        """
        Update the object with the given keyword arguments.

        Parameters
        ----------
        **kwargs : dict
            The keyword arguments to update the object with.

        Returns
        -------
        None
        """
        self.__dict__.update(kwargs)
    
    def Component2Index(self,k,c):
        """
        Returns the index of the component c for variable k in the data array.

        Parameters
        ----------
        k : str
            Variable name.
        c : str
            Component name.

        Returns
        -------
        int
            Index of the component c for variable k in the data array.
        """
        v = self.variables_available[k]
        Component2Index = [v[i] for i in self.components[k]]
        return Component2Index.index(c)
    
    def ParseData(self,data,c):
        """
        Add a capteur to the SEM3DMonitor object.

        Parameters
        ----------
        data : numpy array
            Data to add to the capteur.
        c : int
            Capteur number.

        Returns
        -------
        None
        """
        for k,v in self.variables.items():
            if k != 'Time':
                self.data[k][:,:,c] = data[:,[v[i]+1 for i in self.components[k]*any(self.components[k])+[0]*(not any(self.components[k]))]]
                
    def ParsePosition(self,data,c):
        """
        Add position data for a specific capteur.

        Parameters
        ----------
        data : numpy array
            Array containing position data [x, y, z].
        c : int
            Capteur number.

        Returns
        -------
        None
        """
        self.SetCoordinates(c,data)
            
    def SetCoordinates(self,c,data):
        """
        Set the coordinates for a capteur.

        Parameters
        ----------
        c : int
            Capteur number.
        x, y, z : float
            Coordinates of the capteur.

        Returns
        -------
        None
        """
        self.pos[str(c)] = data[...]
        
    def SetDimensions(self):
        """
        Initialize the dimensions for data and position attributes.

        This method sets up the `data` dictionary with empty numpy arrays for each
        variable (excluding 'Time') based on the specified number of time steps (`nt`),
        components, and capteurs (`nc`). It also initializes the `position` dictionary 
        with empty numpy arrays representing 3D coordinates for each capteur.

        Attributes
        ----------
        data : dict
            A dictionary where each key is a variable name and each value is an
            empty numpy array with dimensions (nt, number of components, nc).
        position : dict
            A dictionary where each key is a capteur index as a string, and each
            value is an empty numpy array representing 3D coordinates.
        """
        for k in self.variables.keys():
            if k!='Time':
                self.data[k] = np.empty((self.nt,max(1,len(self.components[k])),self.nc),
                                        dtype=np.float64)
        # for c in range(self.nc):
        #     self.position[str(c)] = np.empty((3,),dtype=np.float64)         
                
    def SetTime(self,Time):
        """
        Set the time vector `Time` for the SEM3DMonitor object.

        If the format is 'h5' and 'Time' is a variable, sets `nt` to the size of the
        `Time` array. If the time vector is not parsed, raises an exception.
        """
        if not self.TimeSet:
            self.TimeSet = True
            self.Time = Time
            self.setDTime()
            self.SetNt()
        
    def SetNt(self):
        """
        Set the number of time steps `nt` for the SEM3DMonitor object.
        
        If the format is 'h5' and 'Time' is a variable, sets `nt` to the size of the
        `Time` array. If the time vector is not parsed, raises an exception.
        """
        if self.format=='h5' and 'Time' in self.variables:
            try:
                self.nt = self.Time.size
            except:
                raise('Time vector not parsed!')
            
    def setDTime(self):
        """
        Set the time step `dTime` for the SEM3DMonitor object.

        If the format is 'h5' and 'Time' is a variable, calculates `dTime` as
        the difference between the last two elements of the `Time` array. If
        the time vector is not parsed, raises an exception. Otherwise, calculates
        `dTime` based on the data array by taking the midpoint values, accounting
        for potential T0 shifts due to filters.

        Returns
        -------
        float
            The calculated `dTime` value when not using 'h5' format.
        """
        if self.format=='h5' and 'Time' in self.variables.keys():
            try:
                self.dTime = self.Time[-1]-self.Time[-2]
            except:
                raise('Time vector not parsed!')
            
        else:
            # on ne prend pas [1]-[0] car certains filtres decalent le T0
            N=self.data.shape[0]/2
            return self.data[N+1,0] - self.data[N,0]
    
    def CheckComponents(self,v,components):
        return [components_dict[v][c] for c in components if c in components_dict[v].keys() and components_dict[v][c] in self.components[v]]
    def ComputeTotalStress(self, components=['xx','yy','zz','xy','xz','yz']):
        """
        Compute total stresses from deviatoric stresses and pressure.

        Parameters
        ----------
        components : list of str
            Components of the deviatoric stress tensor to use for the computation.
            Default is ['xx','yy','zz','xy','xz','yz'].

        Returns
        -------
        None

        Notes
        -----
        This function assumes that the SEM3DMonitor object has the variables 'StressDev' and
        'Pressure' in its `variables` attribute. If the `components` parameter is not empty,
        the function will create a new variable 'StressTotal' in the `data` attribute of the
        SEM3DMonitor object. The entries of the 'StressTotal' variable are computed as the sum
        of the pressure and the corresponding component of the deviatoric stress tensor.
        """
        if "StressDev" in self.variables.keys() and "Pressure" in self.variables.keys():
            if any(self.components['StressDev']):
                self.data['StressTotal'] = np.empty((self.nt,len(self.CheckComponents("StressDev",components)),self.nc),dtype=np.float64)
                print(f'Computing total stresses for {self.name}...')
                for m in tqdm(range(self.nc)):
                    for c in self.CheckComponents("StressDev",components):
                        if c<3:
                            self.data['StressTotal'] = self.data['Pressure'][:,0,m]+\
                                self.data['StressDev'][:,self.Component2Index('StressDev',c),m]
                        else:
                            self.data['StressTotal'] = self.data['StressDev'][:,self.Component2Index('StressDev',c),m]


    def Plot(self,wkdir,variables,components,monitors,hfg=None,svf=False,**kwargs):
        """
        Plot the time-histories for a given set of capteurs, variables and components.

        Parameters
        ----------
        wkdir : str
            Working directory where the figures are saved if `svf` is True.
        variables : list of str
            List of variables to plot.
        components : list of str
            List of components to plot for each variable.
        monitors : list of int
            List of capteurs to plot.
        hfg : matplotlib.figure.Figure
            If provided, plot on this figure.
        svf : bool
            If True, save the figures in the working directory.
        **kwargs : dict
            Additional arguments to be passed to the plt.plot method.

        Returns
        -------
        fig : matplotlib.figure.Figure
            If `svf` is False, return the figure object.
        """
        if -1 in monitors:
            monitors = range(self.nc)
        print(f'Plotting {self.name}...')
        for m in tqdm(monitors):
            for v in variables:
                if v in self.variables.keys():
                    if any(self.components[v]):
                        for c in self.CheckComponents(v,components):
                            if hfg:
                                plt.figure(hfg.number)
                            else:
                                plt.figure(figsize=[10,5])
                            plt.plot(self.Time,self.data[v][:,self.Component2Index(v,c),m])
                            plt.xlim(self.Time[0],self.Time[-1])
                            plt.xlabel(r'$\mathbf{t [s]}$',fontsize=14)
                            plt.ylabel(r'$\mathbf{{ {{{vv}}} [{{{uu}}}] }}$'.format(vv=v,uu=units[v]),fontsize=14)
                            if svf:
                                plt.savefig(osj(wkdir,self.name+'_'+str(m)+'_'+str(v)+'_'+str(c)+'.png'),\
                                            dpi=300,bbox_inches='tight')
                                plt.close()
                            else:
                                return plt.gcf()
                    else:
                        if hfg:
                            plt.figure(hfg.number)
                        else:
                            plt.figure(figsize=[10,5])
                        plt.plot(self.Time,self.data[v][:,0,m])
                        plt.xlim(self.Time[0],self.Time[-1])
                        plt.xlabel(r'$\mathbf{t [s]}$',fontsize=14)
                        plt.ylabel(r'$\mathbf{{ {{{vv}}} [{{{uu}}}] }}$'.format(vv=v,uu=units[v]),fontsize=14)
                        if svf:
                            plt.savefig(osj(wkdir,self.name+'_'+str(m)+'_'+str(v)+'.png'),\
                                        dpi=300,bbox_inches='tight')
                            plt.close()
                        else:
                            return plt.gcf()
                        
                        # raise ValueError('Variable '+v+' not parsed!')

def ParseSEM3DH5Traces(wkdir='./',
                       format='h5',
                       variables=[''],
                       components=['x','y','z'],
                       names='all',**kwargs):    
    """
    Parse SEM3D HDF5 traces.

    Parameters
    ----------
    wkdir : string
        Working directory where files are located.
    format : string
        Format of traces. Only 'h5' is supported.
    variables : list of strings
        List of variables to parse. If empty, all variables are parsed.
    components : list of strings
        Components to parse. If empty, all components are parsed.
    names : string or list of strings
        Name of the capteur to parse. If 'all', all capteurs are parsed.

    Returns
    -------
    SEM3Dstream : dictionary
        A dictionary of SEM3DMonitor objects, one per capteur.
    """
    # handle h5 format
    if format == 'h5': 
        # get file names
        FileList = glob.glob(osj(wkdir,'*.h5'))
        # read the first file to get variable names and available components
        FileName = FileList[0]
        FileHandle = h5py.File(FileName,"r+")
        # create ordered dictionary for available and selected variables
        variables_available = OrderedDict({})
        variables_ok = OrderedDict({'Time':0})
        # create list of tuples (variable, component)
        tmp = [[v.decode("utf-8") for v in a.split()] for a in FileHandle['Variables'][...].tolist()]
        # check available variables
        for vi,v in enumerate(tmp):
            if len(v)>2:
                v=(v[0]+v[1],v[-1])
            if v[0] not in variables_available.keys():
                variables_available[v[0]] = [int(v[-1])-1]
            else:
                variables_available[v[0]].append(int(v[-1])-1)
        # select variables from available ones
        for vi,v in enumerate(tmp):
            if len(v)>2:
                v=(v[0]+v[1],v[-1])
            if v[0] in variables:
                if v[0] in variables_available.keys():
                    if v[0] not in variables_ok.keys():
                        variables_ok[v[0]] = [vi-1]
                    else:
                        variables_ok[v[0]].append(vi-1)
                else:
                    variables.remove(v[0])
                    print(f"Variable {v[0]} is not available")
                    pass
        # define selected time-history components for each variable
        components_ok = {}
        for k in variables_ok.keys():
            if k!='Time':
                components_ok[k] = [components_dict[k][c] for c in components if c in components_dict[k]]
                components_ok[k].sort()

        # define selected monitor subsets to parse
        all_names = []
        for FileName in FileList:
            FileHandle = h5py.File(FileName,"r+")
            for H5Dataset in FileHandle.items():
                if H5Dataset[0][-4:] == "_pos":
                    if H5Dataset[0][:-4].split('_')[-2] not in all_names:
                        all_names.append(H5Dataset[0][:-4].split('_')[-2])
            FileHandle.close()
        if 'all' in names: 
            names_ok = [c for c in all_names]
        else:
            names_ok = names
        SEM3Dstream = {}
        # Parse h5 files to find total number of capteurs

        for n in names_ok:
            # number of monitors per set (named n)
            nc = 0
            # dictionary pointing each monitor to the belonging file
            idx = {}
            monitor = SEM3DMonitor(name=n, format=format,
                                   variables=variables_ok, 
                                   variables_available=variables_available,
                                   components=components_ok)
            for FileName in FileList:
                idx[FileName] = []
                FileHandle = h5py.File(FileName,"r+")
                for H5Dataset in FileHandle.items():
                    if 'Variables' not in H5Dataset[0]:
                        if "_pos" not in H5Dataset[0]:
                            if any(set(H5Dataset[0].split('_')[:-1]) & set(names_ok)):
                                monitor.SetTime(H5Dataset[1][...][:,0])
                        else:
                            if any(set(H5Dataset[0].split('_')[:-1]) & set(names_ok)):
                                nc += 1 # increase number of found monitors
                                # assign monitor position and belonging filename
                                if len(H5Dataset[0].split('_'))>2:
                                    idx[FileName].append(int(H5Dataset[0].split('_')[-2]))
                                    monitor.ParsePosition(H5Dataset[1],int(H5Dataset[0].split('_')[-2]))
                                else:
                                    idx[FileName].append(0)
                                    monitor.ParsePosition(H5Dataset[1],0)
                FileHandle.close()
            # set index and number of monitoring points
            monitor(nc=nc,idx=idx)
            # set the data storage dimensions
            monitor.SetDimensions()
            # read datasets from h5 files
            for FileName in FileList:
                FileHandle = h5py.File(FileName,"r+")
                for H5Dataset in FileHandle.items():
                    if 'Variables' not in H5Dataset[0] and "_pos" not in H5Dataset[0]:
                        if any(set(H5Dataset[0].split('_')[:-1]) & set(names_ok)):
                            if len(H5Dataset[0].split('_'))>1:
                                monitor.ParseData(H5Dataset[1],int(H5Dataset[0].split('_')[-1]))
                            else:
                                monitor.ParseData(H5Dataset[1],0)
                FileHandle.close()
            SEM3Dstream[n]=monitor
        return SEM3Dstream
    else:
        raise(f"ERROR: format {format} not implemented!")

def ParseOptions():
    """
    Parse command line arguments for ParseSEM3DH5Traces.

    Returns
    -------
    dict
        A dictionary containing the parsed command line arguments.
    """
    OptionParser = argparse.ArgumentParser(prefix_chars='@')
    OptionParser.add_argument('@w','@@wkdir',type=str,default='../test/traces/',help='Database main directory')
    OptionParser.add_argument('@f','@@format',type=str,default='h5',help='Database format [h5|txt]')
    OptionParser.add_argument('@n','@@names',type=str,nargs='+',default=['all'],help = 'Name of the set(s) of monitors')
    OptionParser.add_argument('@v','@@variables',type=str,nargs='+',default=['Displ','Veloc','Accel'],help='Output variables [Displ|Veloc|Accel|Pressure|StressDev|EpsDev|EpsDevPl|EpsDevSv]')
    OptionParser.add_argument('@c','@@components',type=str,nargs='+',default=['x','y','z'],help='Motion components [x|y|z]')
    OptionParser.add_argument('@m','@@monitors',type=int,nargs='+',default=[-1],help='Monitor number. -1 plots all monitors')
    OptionParser.add_argument('@s','@@stress',action='store_true',default=False,help='Compute Total stress?')
    OptionParser.add_argument('@p','@@plot',action='store_true',default=False,help='Plot?')
    if len(sys.argv)==1:
        OptionParser.print_help(sys.stderr)
        sys.exit(1)

    options = OptionParser.parse_args().__dict__
    return options

def main():
    options = ParseOptions()
    
    print("Parse {} database ({}*.{})\nVariables: {}-Comp: {}".format(options['names'],options['wkdir'],options['format'],
                                                              options['variables'],options['components']))
    
    # Parse SEM3D database of simulated time-histories
    stream = ParseSEM3DH5Traces(**options)
    print("Database parsed!\n")

    # Compute total stresses if present
    if options['stress']:
        if any(['StressDev' in list(stream.values())[0].variables_available.keys() and 'Pressure' in list(stream.values())[0].variables_available.keys() for n in options['names']]):
            for n,st in stream.items():
                st.ComputeTotalStress(components=options['components'])
    else:
        raise("Warning: no stress field defined in database")

    # Plot traces
    if options['plot']:
        if 'all' in options['names']:
            for n,st in stream.items():
                st.Plot(**options,svf=True)
        else:
            for n in options['names']:
                stream[n].Plot(**options,svf=True)


if __name__=='__main__':
    main()    
    
