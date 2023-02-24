# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to parse the USGS finite fault files (.fsp) and extract the 
slip patch in a csv file (per segment)


Ex.1 : Parse and plot the generic complete_inversion.fsp file from USGS website
        
        python3 usgs2rik.py @p
"""
import argparse
import pdb
import os
import numpy as np
import dask.dataframe as dd 
from itertools import groupby
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt

# General information
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2023, LMPS UMR CNRS 9026 - CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

class USGSFiniteFault(object):
    def __init__(self, **kwargs):
        super(USGSFiniteFault, self).__init__()
        self.__dict__.update(**kwargs)
        self.options = {}
        self.parsed = False
        self.__fault = {}
        self.read_da = False
        self.ParseCL()
        self.__call__()
        
    def __call__(self):
        self.fault={}
        
    def ParseCL(self): 
        parser = argparse.ArgumentParser(prefix_chars='@')
        parser.add_argument("@w","@@workingdir",
                            type=str,
                            default=r"./",
                            help='Workin directory')
        parser.add_argument("@f","@@filename", 
                            type=str,
                            default=r"complete_inversion.fsp",
                            help='Complete inversion filename')
        parser.add_argument("@n","@name",
                            type=str,
                            default=r"fault_name",
                            help='Fault name')
        parser.add_argument("@p", "@@plot",
                            action='store_true',
                            default=False,
                            help='Plot fault segments?')
        self.options = parser.parse_args().__dict__
        self.options["fullpath"] = os.path.join(self.options["workingdir"],
                                                "{:>s}".format(self.options["filename"])
                                                )

        self.parsed = True
        self.__dict__.update(self.options)        

    def get_comment(self,ln):
        return "%" in ln

    def get_segment(self,ln):
        return "% SEGMENT #" in ln

    def get_lenwid(self,ln):
        return "% LEN =" in ln

    def get_z2top(self,ln):
        return "% depth to top: Z2top =" in ln

    def get_hypo_ll(self,ln):
        return "% LAT =" in ln

    def get_hypo_along(self,ln):
        return "% hypocenter on SEG #" in ln

    def get_nsbfs(self,ln):
        return "% Nsbfs =" in ln

    def get_NxNz(self,ln):
        return "% Invs : Nx =" in ln

    def get_DxDz(self,ln):
        return "% Invs : Dx =" in ln
    
    def SetFiniteFault(self, *args):
        with open(self.fullpath, 'r') as fobj:
            segments = [[ln.split() for ln in group]
                        for has_data, group in groupby(fobj, self.get_segment)
                        if has_data]
            self.__fault['ns'] = len(segments)
            self.__fault["tag"] = np.empty((self.__fault['ns'],), dtype=np.int64)
            self.__fault["stk"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["dip"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["len"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["wid"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["z2t"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["hla"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["hlo"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["hal"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["hdd"] = np.empty((self.__fault['ns'],), dtype=np.float64)
            self.__fault["nsf"] = np.empty((self.__fault['ns'],), dtype=np.int64)
            self.__fault["nx"] = np.empty((1,), dtype=np.int64)
            self.__fault["nz"] = np.empty((1,), dtype=np.int64)
            self.__fault["dx"] = np.empty((1,), dtype=np.float64)
            self.__fault["dz"] = np.empty((1,), dtype=np.float64)

            for i,s in enumerate(segments):
                self.__fault['tag'][i] = float(s[0][3].strip(":"))
                self.__fault['stk'][i] = float(s[0][6])
                self.__fault['dip'][i] = float(s[0][10])
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            lenwid = [[ln.split() for ln in group]
                      for has_data, group in groupby(fobj, self.get_lenwid)
                      if has_data]
            for i, s in enumerate(lenwid):
                self.__fault["len"][i] = float(s[0][3])*1000.0
                self.__fault["wid"][i] = float(s[0][7])*1000.0
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            z2top = [[ln.split() for ln in group]
                     for has_data, group in groupby(fobj, self.get_z2top)
                     if has_data]
            for i, s in enumerate(z2top):
                self.__fault["z2t"][i] = float(s[0][6])*1000.0
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            hypo_ll = [[ln.split() for ln in group]
                       for has_data, group in groupby(fobj, self.get_hypo_ll)
                       if has_data]
            for i, s in enumerate(hypo_ll):
                self.__fault["hla"][i] = float(s[0][3].strip(","))
                self.__fault["hlo"][i] = float(s[0][6].strip(","))
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            hypo_al = [[ln.split() for ln in group]
                       for has_data, group in groupby(fobj, self.get_hypo_along)
                       if has_data]
            for i, s in enumerate(hypo_al):
                self.__fault["hal"][i] = float(s[0][10].strip(","))
                self.__fault["hdd"][i] = float(s[0][14].strip(","))
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            nsbfs = [[ln.split() for ln in group]
                     for has_data, group in groupby(fobj, self.get_nsbfs)
                     if has_data]
            for i, s in enumerate(nsbfs):
                self.__fault["nsf"][i] = int(s[0][3])
            fobj.close()

        with open(self.fullpath, 'r') as fobj:
            DxDz = [[ln.split() for ln in group]
                    for has_data, group in groupby(fobj, self.get_DxDz)
                    if has_data]
            for i, s in enumerate(DxDz):
                self.__fault["dx"][i] = float(s[0][5])*1000.0
                self.__fault["dz"][i] = float(s[0][9])*1000.0
            fobj.close()

        self.__fault["nz"] = (self.__fault["wid"]/self.__fault["dz"]).astype(np.int64)
        self.__fault["nx"] = (self.__fault["len"]/self.__fault["dx"]).astype(np.int64)

        
    def GetFiniteFault(self):
        return self.__fault
    
    fault = property(GetFiniteFault, SetFiniteFault)

    def ParseSlipPatches(self):
        df = dd.read_csv(self.fullpath,
                         header=None,
                         names=["LAT",
                                "LON",
                                "EW",
                                "NS",
                                "Z",
                                "SLIP",
                                "RAKE",
                                "TRUP",
                                "RISE",
                                "SF_MOMENT"],
                         skiprows=58,
                         sep="\s+|;|:",
                         engine="python",
                         comment="%",
                         )
        df.divisions = tuple([0]+[self.fault["nsf"].sum().tolist()])
        df = df.repartition(divisions=tuple([0]+self.fault["nsf"].cumsum().tolist()))
        self.chunk_sizes = list(df.map_partitions(len).compute().values)
        self.slip_da = df[["EW", "NS", "Z", "SLIP"]].to_dask_array(lengths=self.chunk_sizes)
        self.read_da = True
        
    def GenerateRIKfiles(self):
        if not self.read_da:
            self.ParseSlipPatches()
            self.read_da = True
            
        fs = [b.compute() for b in self.slip_da.blocks]
        fs = [np.flip(fs[i][:, 3].reshape((self.fault["nz"][i], 
                                           self.fault["nx"][i])), 
                       axis=0)
               for i in range(len(fs))]

        for i, s in enumerate(fs):
            np.savetxt(os.path.joint(self.fullpath,
                                     "{:s}_segment{:d}.csv".format(self.name,i)
                                    ), 
                       s, 
                       delimiter="\t", 
                       fmt="%7.3f")
        
        np.savetxt(os.path.joint(self.fullpath,
                                 "{:s}_segment_all.csv".format(self.name)
                                ),
                   np.hstack(fs),
                   delimiter="\t",
                   fmt="%7.3f"
                   )
    def PlotFault(self):
        if self.plot:
            fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

            fs = [b.compute() for b in self.slip_da.blocks]

            ls = LightSource(270, 45)
            
            m = cm.ScalarMappable(cmap=cm.jet)
            for i in range(len(fs)):
                rgb = ls.shade(fs[i][:, 3].reshape((self.fault["nz"][i],
                                                    self.fault["nx"][i])),
                               cmap=cm.jet,
                               vert_exag=0.1,
                               blend_mode='soft')
                surf = ax.plot_surface(fs[i][:, 0].reshape((self.fault["nz"][i],
                                                            self.fault["nx"][i])),
                                       fs[i][:, 1].reshape((self.fault["nz"][i],
                                                            self.fault["nx"][i])),
                                       np.flip(fs[i][:, 2].reshape((self.fault["nz"][i],
                                                                    self.fault["nx"][i])),
                                               axis=0),
                                       rstride=1,
                                       cstride=1,
                                       facecolors=rgb,
                                       linewidth=0,
                                       antialiased=False,
                                       shade=False)
            # cbar = fig.colorbar(surf)
            ax.set(xlabel=r"Lon [$^\circ$]", 
                   ylabel=r"Lat [$^\circ$]",
                   zlabel=r"z [m]")
            plt.savefig(os.path.join(self.workingdir,
                                     "{:s}.png".format(self.name)), 
                        bbox_inches='tight')
            plt.close()

if __name__=="__main__":
    
    Fault = USGSFiniteFault(name="fault_name")
    
    Fault.GenerateRIKfiles()
    
    Fault.PlotFault()
