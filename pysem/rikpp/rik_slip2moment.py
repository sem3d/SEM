# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to convert RIKsrf2 output into SEM3D input

Example:

    python3 rik_slip2moment.py @nL 140 @nW 80 @sf ./slipdistribution.dat @hL 3.66 @hW 3.15 @L 7.0 @W 4.0 @fg teil @strike 45. @dip 60. @rake 108. @nt 2048 @dt 0.01 @hE 631892.2 @hN 4931475.0 @hZ -860.0 @wkd ./ @tag teil


"""

#=======================================================================
# Required modules
#=======================================================================
from rik_pp_lib import *
from scipy.integrate import cumtrapz
from scipy.spatial import Voronoi,Delaunay,voronoi_plot_2d
u'''General informations'''
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, MSSMat UMR CNRS 8579 - CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
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
        surf = ax.plot_surface(self.xg, self.yg, self.zg, \
                               rstride=1, cstride=1, cmap=cm.coolwarm,\
                               linewidth=0, antialiased=False)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
    
    def genmesh2d(self):
        p2d = np.array([self.xg.T.reshape(-1,),\
                        self.yg.T.reshape(-1,)],dtype=np.float64).T
        self.msh = Delaunay(p2d)
        self.msh.nodes = np.zeros((self.msh.points.shape[0],self.msh.points.shape[1]+1))
        self.msh.nodes[:,:-1]=self.msh.points
        return 
    
    def genmesh3d(self):
        p2d = np.array([self.xg.T.reshape(-1,),\
                        self.yg.T.reshape(-1,)],dtype=np.float64).T
        zg = self.zg.T.reshape(-1,)
        self.msh = Delaunay(p2d)
        self.msh.nodes = np.zeros((self.msh.points.shape[0],\
                                   self.msh.points.shape[1]+1),\
                                   dtype=np.float64)
        self.msh.nodes[:,:-1]=self.msh.points 
        for s in self.msh.simplices:
            self.msh.nodes[s,-1] = zg[s]
        return
    
    def RotTransMesh2d(self,aS,aD,hyp):
        MatMesh = get_rotation_tensor(aS, aD)
        trs = np.dot(MatMesh,hyp[0])
        trs = hyp[1]-trs
        #self.msh.nodes = np.empty((self.msh.points.shape[0],self.msh.points.shape[1]+1))
        for i,p in enumerate(self.msh.points):
            self.msh.nodes[i,:] = trs+np.dot(MatMesh,\
                                             np.concatenate((p.flatten(),\
                                                             np.array([0.],dtype=np.float64))))
            
    def RotTransMesh3d(self,aS,aD,hyp):
        MatMesh = get_rotation_tensor(aS, aD)
        trs = np.dot(MatMesh,hyp[0])
        trs = hyp[1]-trs
        for i,p in enumerate(self.msh.nodes):
            self.msh.nodes[i,:] = trs+np.dot(MatMesh,p.flatten())
            
    def write_mesh2h5(self,fid):
        fid.create_dataset(name='Nodes', data=self.msh.nodes)
        fid.create_dataset(name='Elements', data=self.msh.simplices)

class SEM_source(object):
    def __init__(self,e=0.,n=0.,z=0.):
        self.e = e
        self.n = n
        self.z = z
    def __call__(self,e=0.,n=0.,z=0.):
        self.e = e
        self.n = n
        self.z = z
    def get_hypo(self):
        return np.array([self.e,self.n,self.z],dtype=np.float64)

class RIK_source(object):
    def __init__(self,x=0.,y=0.,z=0.,hypfile=None):
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
        return np.array([self.x,self.y,self.z],dtype=np.float64)

if __name__=='__main__':
    opt = start_rik()
    globals().update(opt)
    fg = opj(wkd,"{:>s}".format(tag))
    coord_file_name = 'source_coordinates.csv'

    # Input files 
    RIK_slipfile = opj(wkd,'slipdistribution.dat')
    mrfile = opj(wkd,'MomentRate.dat')
    srfile = opj(wkd, 'sr.dat')
    hypfile = opj(wkd,'nucleationpoint.dat')

    T_total = (nt-1)*dt
    time = np.linspace(0.0, T_total, nt)
    # Les angles en radian
    print("/n/n")
    print("Hypocenter (UTM): \nEW: {:>.1f}\nNS: {:>.1f}\nUD: {:>.1f}\n".format(hE,hN,hZ))
    print("strike: {:>.1f}, dip: {:>.1f}, rake: {:>.1f}\n".format(strike,dip,rake))
    aS,aR,aD = np.pi/180.0*strike,np.pi/180.0*rake,np.pi/180.0*dip

    MatMesh = get_rotation_tensor(aS,aD)
    src_sem = SEM_source(e=hE,n=hN,z=hZ)
    src_rik = RIK_source(hypfile=hypfile,z=np.abs(hZ)/1.0e3)
    RIK_hypocenter,SEM_hypocenter = src_rik.get_hypo(),src_sem.get_hypo()
    print("RIK_hypocenter: \nx [km]: {0}\ny [km]: {1}\nd [km]: {2}\n".format(*tuple(RIK_hypocenter.tolist())))
    MR = read_moment_rate_RIK(mrfile,(int(nL),int(nW),int(nt)))
    SR = read_moment_rate_RIK(srfile,(int(nL),int(nW),int(nt)))
    
    # Creation des fichiers d'entree pour SEM
    kfo = h5py.File(opj(wkd,'{:>s}_kine.h5'.format(tag)),'w')
    sfo = h5py.File(opj(wkd,'{:>s}_moment.h5'.format(tag)),'w')
    #driver='mpio', comm=MPI.COMM_WORLD)
    # Attribuer les proprietes temporelles
    kfo.attrs['Nt'] = nt
    kfo.attrs['dt'] = time[1]-time[0]
    kfo.attrs['Ns'] = nL 
    kfo.attrs['Nd'] = nW
    sfo.create_dataset('time', data=time)
    # Vecteurs du normale et du slip
    vecteurs(aD,aS,aR,kfo)

    # Coordonnees des points dans le repere SEM3D
    xgrid = kfo.create_dataset('x',(nL,nW),chunks=(1,1))
    ygrid = kfo.create_dataset('y',(nL,nW),chunks=(1,1))
    depth = kfo.create_dataset('z',(nL,nW),chunks=(1,1))
    xgrid[:,:] = 0.0
    ygrid[:,:] = 0.0
    depth[:,:] = 0.0

    RIK_xcoord = np.genfromtxt(RIK_slipfile, usecols=1)  # Y_RIK
    RIK_ycoord = np.genfromtxt(RIK_slipfile, usecols=0)  # X_RIK

    # Creer le matrice hdf5 pour le moment
    Moment = sfo.create_dataset('moment',(nL,nW,nt),chunks=(1,1,nt))
    Slip=np.zeros((nL,nW,nt),dtype=np.float32)
    # Integration pour calculer le moment
    n = 0
    for i in range(0, nL):
        for j in range(0, nW):
            n = n+ 1
            Moment[i,j,:] = cumtrapz(MR[i,j,:],dx=dt,initial=0.)
            Slip[i,j,:] = cumtrapz(SR[i,j,:],dx=dt,initial=0.)
            
    for t in range(nt):
        sfo.create_dataset('mom_{:>d}'.format(t), shape=(nL*nW,),\
                           data = Moment[:,:,t].T.reshape((nL*nW,)),\
                           chunks=(1,))
        kfo.create_dataset('slp_{:>d}'.format(t), shape=(nL*nW,),\
                           data = Slip[:,:,t].T.reshape((nL*nW,)),\
                           chunks=(1,))
        kfo.create_dataset('sra_{:>d}'.format(t), shape=(nL*nW,),\
                           data = SR[:,:,t].T.reshape((nL*nW,)),\
                           chunks=(1,))
    
    n = 0
    for j in np.arange(nW):
        for i in np.arange(nL):
            if (nL*nW > 1):
                xgrid[i,j] = RIK_xcoord[n]
                ygrid[i,j] = RIK_ycoord[n]
            else:
                xgrid[i,j] = RIK_xcoord
                ygrid[i,j] = RIK_ycoord
            depth[i,j] = RIK_hypocenter[2]+np.sin(aD)*(RIK_hypocenter[1]-xgrid[i,j])
            xgrid[i,j] *= 1.0e3 
            ygrid[i,j] *= 1.0e3 
            depth[i,j] *= 1.0e3
            n = n+1
    RIK_hypocenter*=1.0e3
    fmsh = mesh(xg=xgrid[()],yg=ygrid[()],zg=depth[()])
    fmsh.RotTransMesh3d(aS,aD,(RIK_hypocenter,SEM_hypocenter))

    fmsh.write_mesh2h5(sfo)
    fmsh.write_mesh2h5(kfo)
    exit()
    fark  = np.dot(MatMesh,RIK_hypocenter)
    fark  = (SEM_hypocenter-fark)
    # Opening file to save SEM3D coordinates of grid points
    coord_file = open (coord_file_name, 'w+')
    # Burada rotation yapiyorum
    n = 0
    dx_plane = xgrid[0,1]-xgrid[0,0]
    dy_plane = ygrid[1,0]-ygrid[0,0]
    dz_plane = depth[0,1]-depth[0,0]
    dx_rot = np.dot(MatMesh,np.array([dx_plane,dy_plane,dz_plane]))
    coord_file.write('{:>10},{:>15},{:>15},{:>15}\n'.format('N','X','Y','Z'))
    for j in np.arange(nW):
        for i in np.arange(nL):
            n = n+ 1
            # Rotasyon uyguluyorum
            coord = np.array([xgrid[i,j], ygrid[i,j], depth[i,j]])
            dum   = np.dot(MatMesh,coord)

            # Aradaki farki ekliyorum; boylelikle translate ediyor
            xgrid[i,j] = dum[0]+ fark[0]
            ygrid[i,j] = dum[1]+ fark[1]
            depth[i,j] = dum[2]+ fark[2]

            # Writing out the cordinates
            coord_file.write('{:>10d},{:>15.5f},{:>15.5f},{:>15.5f}\n'.format(\
                n,xgrid[i,j],ygrid[i,j],depth[i,j]))
    coord_file.close()

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
        ax.set_xlabel(r'$t$ [s]'   ,fontsize=20)
        ax.set_ylabel(r'$M_0$ [Nm]',fontsize=20)
        plt.gca().spines["top"].set_alpha(0)
        plt.gca().spines["bottom"].set_alpha(.3)
        plt.gca().spines["right"].set_alpha(0)
        plt.gca().spines["left"].set_alpha(.3)
        # ax.legend()
        fig.savefig("{:>s}_moment_ths.png".format(fg),dpi=500,bbox_inches='tight')
        fig.savefig("{:>s}_moment_ths.eps".format(fg),dpi=500,bbox_inches='tight',
                format='eps')
        plt.close()
    # Fermeture des fichiers hdf5
    kfo.close()
    sfo.close()
