# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to convert RIKsrf2 output into SEM3D input

Example:

    python3 rik_slip2moment.py @nL 140 @nW 80 @sf ./slipdistribution.dat 
    @hL 3.66 @hW 3.15 @L 7.0 @W 4.0 @fg teil @strike 45. @dip 60. @rake 108. 
    @nt 2048 @dt 0.01 @hE 631892.2 @hN 4931475.0 @hZ -860.0 @wkd ./ @tag teil


"""
# Required modules
from rik_pp_lib import *
from scipy.integrate import cumtrapz
from scipy.spatial import Voronoi,Delaunay,voronoi_plot_2d

# General informations
__author__ = "Filippo Gatti and Elif Oral"
__copyright__ = "Copyright 2018, MSSMat UMR CNRS 8579, CentraleSupélec"
__credits__ = ["Filippo Gatti", "Elif Oral"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

class mesh(object):
    def __init__(self,**args):
        if len(args.keys())==1:
            self.xg = args['xg']
        if len(args.keys())==2:
            self.xg = args['xg']
            self.yg = args['yg']
        if len(args.keys())==3:
            self.xg = args['xg']
            self.yg = args['yg']
            self.zg = args['zg']
        if self.xg.size:
            self.nW = self.xg.shape[0]
            self.nL = self.xg.shape[1]
        if self.xg.size and self.yg.size:
            try:
                self.genmesh3d()
            except:
                self.genmesh2d()
               
    def __call__(self,**args):
        self.__init__(**args)
    
    def plot_3d_geo(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(self.xg, 
                               self.yg, 
                               self.zg, \
                               rstride=1, 
                               cstride=1, 
                               cmap=cm.coolwarm,\
                               linewidth=0, 
                               antialiased=False)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
    
    def genmesh2d(self):
        p2d = np.array([self.xg.T.reshape(-1,),\
                        self.yg.T.reshape(-1,)],dtype=np.float64).T
        self.msh = Delaunay(p2d)
        self.msh.nodes = np.zeros((self.msh.points.shape[0],
                                   self.msh.points.shape[1]+1))
        self.msh.nodes[:,:-1]=self.msh.points
        return 
    
    def genmesh3d(self):
        p2d = np.array([self.xg.T.reshape(-1,),\
                        self.yg.T.reshape(-1,)],
                       dtype=np.float64).T
        zg = self.zg.T.reshape(-1,)
        self.msh = Delaunay(p2d)
        self.msh.nodes = np.zeros((self.msh.points.shape[0],\
                                   self.msh.points.shape[1]+1),\
                                   dtype=np.float64)
        self.msh.nodes[:,:-1]=self.msh.points 
        for s in self.msh.simplices:
            self.msh.nodes[s,-1] = zg[s]
        return
    
    def RotTransMesh2d(self,φs,δ,hyp):
        Q_φδλ = get_rotation_tensor(φs, δ)
        trs = np.dot(Q_φδλ,hyp[0])
        trs = hyp[1]-trs
        #self.msh.nodes = np.empty((self.msh.points.shape[0],self.msh.points.shape[1]+1))
        for i,p in enumerate(self.msh.points):
            self.msh.nodes[i,:] = trs+np.dot(Q_φδλ,\
                                             np.concatenate((p.flatten(),\
                                                             np.array([0.],
                                                                      dtype=np.float64))))
            
    def RotTransMesh3d(self,φs,δ,hyp):
        Q_φδλ = get_rotation_tensor(φs, δ)
        trs = np.dot(Q_φδλ,hyp[0])
        trs = hyp[1]-trs
        for i,p in enumerate(self.msh.nodes):
            self.msh.nodes[i,:] = trs+np.dot(Q_φδλ,
                                             p.flatten())
            
    def write_mesh2h5(self,fid):
        fid.create_dataset(name='Nodes', 
                           data=self.msh.nodes)
        fid.create_dataset(name='Elements', 
                           data=self.msh.simplices)


# class classproperty(object):
#     def __init__(self, fget):
#         self.fget = classmethod(fget)

#     def __get__(self, obj, owner):
#         return self.fget.__get__(None, owner)()

# class MyCalculatedClassAttrsMeta(type):
#     def __new__(cls, name, bases, dct):
#         c = super().__new__(cls, name, bases, dct)
#         c.__nSegments += 1
#         return c


class SEM3Dfault(object):
    __nSegments = 0
    __hypocenter = np.empty((3,))
    __M0 = 0.0
    
    def __init__(self, **kwargs):
        super(SEM3Dfault, self).__init__()
        if 'hypocenter' in kwargs.keys():
            SEM3Dfault.__hypocenter = kwargs['hypocenter']
        if 'moment' in kwargs.keys():
            SEM3Dfault.__M0 = kwargs['moment']
            
    def __call__(self, **kwargs):
        if 'hypocenter' in kwargs.keys():
            SEM3Dfault.__hypocenter = kwargs['hypocenter']
        if 'moment' in kwargs.keys():
            SEM3Dfault.__M0 = kwargs['moment']
            
    @classproperty
    def UpdateNumberOfSegments(cls):
        cls.__nSegments += 1
    
    # def GetNumberOfSegments(self):
    #     return SEM3Dfault.__nSegments
    
    def SetHypoCoords(SEM3Dfault, hypocenter):
        for k in sorted(hypocenter, 
                        key=hypocenter.get, 
                        reverse=True):
            print("{:>10}{:>20}".format(k, hypocenter[k]))
            
        hypocenter2array=np.array([v for v in list(hypocenter.values())])
        SEM3Dfault.__hypocenter = hypocenter2array

    def GetHypoCoords(SEM3Dfault):
        return SEM3Dfault.__hypocenter
    
    HypoCoords = property(GetHypoCoords, SetHypoCoords)
    
    @classproperty
    def M0(cls):
        return cls.__M0
    
    @M0.setter
    def SetM0(cls, M0):
        cls.__M0 = M0
        
    @classproperty
    def nSegments(cls):
        return cls.__nSegments
    
    @nSegments.setter
    def SetnSegments(cls, count):
        cls.__nSegments += count

class FaultSegment(SEM3Dfault):
    def __init__(self, 
                 strike=None, 
                 dip=None,
                 rake=None):
        super(FaultSegment, self).__init__()
        SEM3Dfault().UpdateNumberOfSegments
        self.SegmentID = SEM3Dfault.nSegments
        self.__φs = None
        self.__δ = None
        self.__λ = None
        self.__sdr = np.empty((3,))
        self.__nv = np.empty((3,))
        self.__dv = np.empty((3,))
        self.__Mv = np.empty((3,3))
        self.set = False
        if strike and dip and rake:
            self.__φs = strike*np.pi/180.0
            self.__δ = dip*np.pi/180.0
            self.__λ = rake*np.pi/180.0
            self.__sdr = np.array([strike, dip, rake])*np.pi/180
            self.set = True
            
    @property
    def SegmentId(self):
        return self.__SegmentID
    
    @property
    def StrikeDipRake(self):
        return self.__sdr

    @property
    def NormalVector(self):
        return self.__nv

    @property
    def SlipVector(self):
        return self.__dv

    @property
    def MomentUnitTensor(self):
        return self.__Mm

    @property
    def FaultVectors(self):
        return self.__ndM
    
    @StrikeDipRake.setter
    def SetStrikeDipRake(self, sdr):
        strike, dip, rake = sdr
        print("Strike: {} Dip: {} Rake: {}".format(strike, dip, rake))
        self.__φs = strike*np.pi/180.0
        self.__δ = dip*np.pi/180.0
        self.__λ = rake*np.pi/180.0
        self.__sdr = np.array([strike, dip, rake])*np.pi/180
        self.set = True
        self.__ndM = compute_seismic_moment_vectors(strike=self.__φs,
                                                    dip=self.__δ,
                                                    rake=self.__λ)
        self.__nv, self.__dv, self.__Mm = self.__ndM
    
    @NormalVector.setter
    def SetNormalVector(self, sdr):
        if not self.set:
            self.SetStrikeDipRake(strike, dip, rake)
        ndM = compute_seismic_moment_vectors(strike=self.__φs,
                                             dip=self.__δ,
                                             rake=self.__λ)
        self.__nv = ndM[0]        
    
    @SlipVector.setter
    def SetSlipVector(self, sdr = None):
        if not self.set:
            self.SetStrikeDipRake(strike, dip, rake)
        ndM = compute_seismic_moment_vectors(strike=self.__φs,
                                             dip=self.__δ,
                                             rake=self.__λ)
        self.__dv = ndM[1]
        
    @MomentUnitTensor.setter
    def SetMomentUnitTensor(self, sdr):
        if not self.set:
            self.SetStrikeDipRake(strike, dip, rake)
        ndM = compute_seismic_moment_vectors(strike=self.__φs,
                                             dip=self.__δ,
                                             rake=self.__λ)
        self.__Mm = ndM[2]

    @FaultVectors.setter
    def SetFaultVectors(self, sdr):
        if not self.set:
            self.SetStrikeDipRake(strike, dip, rake)
        self.__ndM = compute_seismic_moment_vectors(strike=self.__φs,
                                                    dip=self.__δ,
                                                    rake=self.__λ)
        nv, dv, Mm = self.__ndM
        self.SetNormalVector(nv)
        self.SetSlipVector(dv)
        self.MomentUnitTensor(Mm)
        
        
class RIK_source(object):
    def __init__(self, x=0., y=0., z=0., hypfile=None):
        self.x = x
        self.y = y
        self.z = z
        if hypfile:
            self.read_hypo(hypfile)
    def __call__(self,x=0.,y=0.,z=0.):
        self.x = x
        self.y = y
        self.z = z
        if hypfile:
            self.read_hypo(hypfile)
    def read_hypo(self,hypfile):
        self.x, self.y = readhypo(hypfile)
    def get_hypo(self):
        return np.array([self.x, self.y, self.z],
                        dtype=np.float64)

if __name__=='__main__':
    opt = start_rik()
    globals().update(opt)
    
    fg = opj(wkd,"{:>s}".format(tag))
    coord_file_name = opj(wkd, 'source_coordinates.csv')

    # Input files 
    RIKslipFilename = opj(wkd,'slipdistribution.dat')
    mrfile = opj(wkd,'MomentRate.dat')
    srfile = opj(wkd, 'sr.dat')
    hypfile = opj(wkd,'nucleationpoint.dat')

    T_total = (nt-1)*dt
    time = np.linspace(0.0, T_total, nt)
    
    assert len(strike)==len(dip) and len(dip)==len(rake)
    
    nSegments = len(strike)
    print("Number of segments: {:d}".format(nSegments))
    
    
    SEMsource = SEM3Dfault(moment=M0)
    SEMsource.HypoCoords = {'ew': hE, 'ns': hN, 'ud': hZ}
    import pdb
    pdb.set_trace()
    SEMsegments = []
    for i, (φs, δ, λ) in enumerate(zip(strike, dip, rake)):
        SEMsegments.append(FaultSegment(strike=φs, dip=δ, rake=λ))
    
    
    
    # Location of the hypocenter
    RIKsource = RIK_source(hypfile=hypfile, z=np.abs(hZ)/1.0e3)
    RIKhypocenter, SEMhypocenter = RIKsource.get_hypo(), SEMsource.get_hypo()
    # print("RIKhypocenter: \nx [km]: {0}\ny [km]: {1}\nd [km]: {2}\n".format(
    #     *tuple(RIKhypocenter.tolist())))
    # print("/n/n")
    # print("Hypocenter (UTM): \nEW: {:>.1f}\nNS: {:>.1f}\nUD: {:>.1f}\n".format(hE,
    #                                                                            hN,
    #                                                                            hZ)
    #       )
    RIK_xcoord = np.genfromtxt(RIKslipFilename, usecols=1)
    RIK_ycoord = np.genfromtxt(RIKslipFilename, usecols=0)

    for (i,s), d, r in zip(enumerate(strike), dip, rake):
        print("Segment {:d}:".format(i+1),
              "strike: {: > .1f},".format(s),
              "dip: {: > .1f},".format(d),
              "rake: {: > .1f}\n".format(r)
              )

    # # convert strike, dip, rake into radians
    # φs, δ, λ = map(lambda xx: [np.pi/180.0*x for x in xx], 
    #                (strike, dip, rake))
    
    # # Compute normal and slip vectors
    # nv, dv, Mm = list(map(compute_seismic_moment_vectors, δ, φs, λ))
    
    # # Get rotation tensor per segment
    # Q_φδλ = list(map(get_rotation_tensor, φs, δ))
    
    #  Parse rate files from RIKsrf
    RIKmoment, RIKslip = map(parseRIKrates, 
                 [mrfile, srfile],
                 [(int(nL), int(nW), int(nt))]*2)
    RIK_moment = cumtrapz(RIKmoment, dx=dt, initial=0., axis=-1)
    RIK_slip = cumtrapz(RIKslip, dx=dt, initial=0., axis=-1)
    
    def get_segment_Lidx(x): return int(x/L*nL)
    # def get_segment_Widx(x): return int(x/W*nW)
    
    nLxSegment = list(map(get_segment_Lidx, segL))
    nWxSegment = [nW]*nSegments
    
    SegmentEdgesL = np.cumsum(np.array([0]+nLxSegment))
    SegmentEdgesW = np.array([0]+nWxSegment)
    
    # Create HDF5 files for SEM
    for nLs, nWs, i in zip(nLxSegment,
                           nWxSegment,
                           range(nSegments)
                           ):
        kfo = h5py.File(opj(wkd,'{:>s}_kine_{:d}.h5'.format(tag,i)),'w')
        sfo = h5py.File(opj(wkd,'{:>s}_moment_{:d}.h5'.format(tag,i)),'w')
        #driver='mpio', comm=MPI.COMM_WORLD)
        
        # Assign general properties
        sfo.create_dataset('time', data=time)
        kfo.attrs['dt'] = time[1]-time[0]
        kfo.attrs['Nt'] = nt
        kfo.attrs['Ns'] = nLs
        kfo.attrs['Nd'] = nWs
        kfo.attrs['Vnormal'] = nv[i]
        kfo.attrs['Vslip'] = dv[i]
        
        Moment = sfo.create_dataset('moment',
                                    (nLs, nWs, nt),
                                    chunks=(1, 1, nt)
                                    )
        Slip = np.zeros((nLs, nWs, nt),
                        dtype=np.float32
                        )
        idx = np.ix_(np.linspace(SegmentEdgesL[i],
                                 SegmentEdgesL[i+1],
                                 nLs,
                                 endpoint=False,
                                 dtype=np.int64),
                     np.linspace(SegmentEdgesW[0],
                                 SegmentEdgesW[1],
                                 nWs,
                                 endpoint=False,
                                 dtype=np.int64),
                     np.linspace(0,
                                 nt,
                                 nt,
                                 endpoint=False,
                                 dtype=np.int64),
                     )
        
        Moment[:, :, :] = RIK_moment[idx]
        Slip[:, :, :] = RIK_slip[idx]
        for t in range(time.size):
            sfo.create_dataset('mom_{:>d}'.format(t), shape=(nLs*nWs,),
                               data=Moment[:, :, t].T.reshape((nLs*nWs,)),
                               chunks=(1,))
            kfo.create_dataset('slp_{:>d}'.format(t), shape=(nLs*nWs,),
                               data=Slip[:, :, t].T.reshape((nLs*nWs,)),
                               chunks=(1,))
            kfo.create_dataset('sra_{:>d}'.format(t), shape=(nLs*nWs,),
                               data=RIKslip[:, :, t].T.reshape((nLs*nWs,)),
                               chunks=(1,))
    
        # SEM grid coordinates
        Xg = kfo.create_dataset('x',
                                (nLs, nWs),
                                chunks=(1, 1))
        Yg = kfo.create_dataset('y',
                                (nLs, nWs),
                                chunks=(1, 1))
        depth = kfo.create_dataset('z',
                                   (nLs, nWs),
                                   chunks=(1, 1))
        Xg[:, :] = 0.0
        Yg[:, :] = 0.0
        depth[:, :] = 0.0
    
    n = 0
    import pdb
    pdb.set_trace()
    for j in np.arange(nWs):
        for i in np.arange(nLs):
            if (nLs*nWs > 1):
                Xg[i,j] = RIK_xcoord[n]
                Yg[i,j] = RIK_ycoord[n]
            else:
                Xg[i,j] = RIK_xcoord
                Yg[i,j] = RIK_ycoord
            depth[i,j] = RIKhypocenter[2]+np.sin(δ)*(RIKhypocenter[1]-Xg[i,j])
            Xg[i,j] *= 1.0e3 
            Yg[i,j] *= 1.0e3 
            depth[i,j] *= 1.0e3
            n = n+1
    RIKhypocenter*=1.0e3
    fmsh = mesh(xg=Xg[()], yg=Yg[()], zg=depth[()])
    fmsh.RotTransMesh3d(φs, δ, (RIKhypocenter, SEMhypocenter))

    fmsh.write_mesh2h5(sfo)
    fmsh.write_mesh2h5(kfo)
    # exit()
    # fark  = np.dot(Q_φδλ,RIKhypocenter)
    # fark  = (SEMhypocenter-fark)
    # # Opening file to save SEM3D coordinates of grid points
    # coord_file = open (coord_file_name, 'w+')
    # # Burada rotation yapiyorum
    # n = 0
    # dx_plane = Xg[0,1]-Xg[0,0]
    # dy_plane = Yg[1,0]-Yg[0,0]
    # dz_plane = depth[0,1]-depth[0,0]
    # dx_rot = np.dot(Q_φδλ,np.array([dx_plane,dy_plane,dz_plane]))
    # coord_file.write('{:>10},{:>15},{:>15},{:>15}\n'.format('N','X','Y','Z'))
    # for j in np.arange(nW):
    #     for i in np.arange(nL):
    #         n = n+ 1
    #         # Rotasyon uyguluyorum
    #         coord = np.array([Xg[i,j], Yg[i,j], depth[i,j]])
    #         dum   = np.dot(Q_φδλ,coord)

    #         # Aradaki farki ekliyorum; boylelikle translate ediyor
    #         Xg[i,j] = dum[0]+ fark[0]
    #         Yg[i,j] = dum[1]+ fark[1]
    #         depth[i,j] = dum[2]+ fark[2]

    #         # Writing out the cordinates
    #         coord_file.write('{:>10d},{:>15.5f},{:>15.5f},{:>15.5f}\n'.format(\
    #             n,Xg[i,j],Yg[i,j],depth[i,j]))
    # coord_file.close()

    if plot:
        clr = sns.color_palette("cool",nL*nW)
        fig = plt.figure(figsize=(10,5))
        sns.set_style('whitegrid')
        ax  = fig.add_subplot(111)
        for i in range(0, nL):
            for j in range(0, nW):
                ax.plot(time,Moment[i,j,:],label='Point '+str(i+1),
                        color=clr[nW*i+j])
        # #
        ax.set_xlabel(r'$t$ [s]',fontsize=20)
        ax.set_ylabel(r'$M_0$ [Nm]',fontsize=20)
        plt.gca().spines["top"].set_alpha(0)
        plt.gca().spines["bottom"].set_alpha(.3)
        plt.gca().spines["right"].set_alpha(0)
        plt.gca().spines["left"].set_alpha(.3)
        # ax.legend()
        fig.savefig("{:>s}_moment_ths.png".format(fg), 
                    dpi=500, 
                    bbox_inches='tight')
        fig.savefig("{:>s}_moment_ths.eps".format(fg), 
                    dpi=500, 
                    bbox_inches='tight',
                format='eps')
        plt.close()
    # Close hdf5
    kfo.close()
    sfo.close()
