from rik_pp_lib import *
from scipy.spatial import Delaunay


class FaultMesh(object):
    def __init__(self, **kwargs) -> None:
        """
        Initialize the class object.

        Parameters
        ----------
        xg : 2D numpy array
            The x coordinates of the grid points on the fault.
        yg : 2D numpy array
            The y coordinates of the grid points on the fault.
        zg : 2D numpy array
            The z coordinates of the grid points on the fault.

        Notes
        -----
        The x and y coordinates are the coordinates of the grid points
        on the fault surface in the local coordinate system of the fault.
        The z coordinates are the depth of the grid points from the surface
        of the Earth.
        """
        self.xg = kwargs['xg']
        self.yg = kwargs['yg']
        self.zg = kwargs['zg']
        if self.xg.size:
            self.nW = self.xg.shape[0]
            self.nL = self.xg.shape[1]
        else:
            ValueError('xg, yg, and zg must be non-empty 2D meshgrids')

    def __call__(self, **kwargs) -> None:
        """
        Initialize the class object.

        This is a convenience method to initialize the class object
        with the same parameters as the __init__ method. This is useful
        when the class object is used as a function.

        Parameters
        ----------
        kwargs : dict
            The keyword arguments to pass to the __init__ method.
        """
        self.__init__(**kwargs)

    def plot_3d_geo(self) -> None:
        """
        Plot the 3D geometry of the fault plane.

        This function will create a 3D plot of the fault plane
        using the x, y, and z coordinates of the grid points.
        The colorbar will represent the depth of the points
        from the Earth's surface.

        Returns
        -------
        None
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(self.xg,
                               self.yg,
                               self.zg,
                               rstride=1,
                               cstride=1,
                               cmap=cm.coolwarm,
                               linewidth=0,
                               antialiased=False)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5,
                     label='Depth (km)')
        plt.show()
    def Triangulate(self) -> None:
        """
        Generate a 3D triangular mesh from the 2D grid points.

        This function will create a 3D triangular mesh from the
        x, y, and z coordinates of the grid points. The mesh
        will be stored in the 'msh' attribute of the class.

        Returns
        -------
        None
        """
        # [CHECK]
        p2d = np.array([self.xg.T.reshape(-1,),
                        self.yg.T.reshape(-1,)],
                       dtype=np.float64).T
        self.msh=Delaunay(p2d)

    def RotTransMesh3d(self, ϕ: float, δ: float, hyp: dict) -> None:
        """
        Rotate and translate the 3D mesh according to the
        fault's orientation and hypocenter location.

        Parameters
        ----------
        ϕ : float
            The strike angle in radians.
        δ : float
            The dip angle in radians.
        hyp : list[float]
            The hypocenter coordinates in the geographical
            coordinate system (East, North, Up).

        Returns
        -------
        None
        """
        # Get rotation tensor
        Q_ϕδλ = get_rotation_tensor(ϕ, δ)
        # Rotate the position around the origin within the fault plane
        self.xg, self.yg, self.zg = np.einsum('ij,jkl->ikl', 
                                              Q_ϕδλ,
                                              np.array([self.xg,
                                                        self.yg,
                                                        self.zg]))
        # Translate the mesh to the hypocenter location
        self.xg += hyp[0]
        self.yg += hyp[1]
        self.zg += hyp[2]
        
        Z = np.zeros((self.msh.points.shape[0],))
        Z = np.row_stack((self.msh.points.T, Z))
        
        self.msh.nodes = np.einsum('ij,jk->ik', Q_ϕδλ, Z).T
        self.msh.nodes +=hyp
        
    def write_mesh2h5(self, fid):
        fid.create_dataset(name='x',
                           data=self.xg)
        fid.create_dataset(name='y',
                           data=self.yg)
        fid.create_dataset(name='z', 
                           data=self.zg)
        fid.create_dataset(name='Nodes',
                           data=self.msh.nodes)
        fid.create_dataset(name='Elements',
                           data=self.msh.simplices)


class SEM3Dfault(object):
    __nSegments = 0
    __HypoXYZ = np.empty((3,))
    __M0 = 0.0
    __SegmentEdgesL = 0.0
    __SegmentEdgesW = 0.0

    def __init__(self, **kwargs):
        super(SEM3Dfault, self).__init__()
        if 'HypoXYZ' in kwargs.keys():
            SEM3Dfault.__HypoXYZ = kwargs['HypoXYZ']
        if 'M0' in kwargs.keys():
            SEM3Dfault.__M0 = kwargs['M0']

    def __call__(self, **kwargs):
        if 'HypoXYZ' in kwargs.keys():
            SEM3Dfault.__HypoXYZ = kwargs['HypoXYZ']
        if 'M0' in kwargs.keys():
            SEM3Dfault.__M0 = kwargs['M0']

    @classproperty
    def UpdateNumberOfSegments(cls):
        cls.__nSegments += 1

    @classproperty
    def HypoXYZ(cls):
        return cls.__HypoXYZ

    @HypoXYZ.setter
    def SetHypoCoords(cls, HypoXYZ: dict):
        for k in sorted(HypoXYZ,
                        key=HypoXYZ.get,
                        reverse=True):
            print("{:>10}{:>20}".format(k, HypoXYZ[k]))

        HypoXYZ = np.array([v for v in list(HypoXYZ.values())])
        cls.__HypoXYZ = HypoXYZ
        print(
            "Hypocenter location: EW: {0:>f} m - NS: {1:>f} m - UD: {2} m\n\n".format(*tuple(cls.HypoXYZ)))

    # HypoXYZ = classproperty(GetHypoCoords, SetHypoCoords)

    @classproperty
    def M0(cls):
        return cls.__M0

    @M0.setter
    def SetM0(cls, M0: float):
        cls.__M0 = M0

    @classproperty
    def nSegments(cls):
        return cls.__nSegments

    @nSegments.setter
    def SetnSegments(cls, count: int):
        cls.__nSegments += count

    @classproperty
    def SegmentEdgesL(cls):
        return cls.__SegmentEdgesL

    @SegmentEdgesL.setter
    def SetSegmentEdgesL(cls, SegmentEdgesL: int) -> None:
        cls.__SegmentEdgesL = SegmentEdgesL

    @classproperty
    def SegmentEdgesW(cls):
        return cls.__SegmentEdgesW

    @SegmentEdgesW.setter
    def SetSegmentEdgesW(cls, SegmentEdgesW: int) -> None:
        cls.__SegmentEdgesW = SegmentEdgesW


def get_segment_Lidx(x, L: float, nL: int): return int(x/L*nL)
# def get_segment_Widx(x): return int(x/W*nW)

class FaultSegment(SEM3Dfault):
    def __init__(self, 
                 strike=None, 
                 dip=None,
                 rake=None):
        super(FaultSegment, self).__init__()
        SEM3Dfault().UpdateNumberOfSegments
        self.SegmentID = SEM3Dfault.nSegments
        self.__ϕ = None
        self.__δ = None
        self.__λ = None
        self.__nLs = 0
        self.__nWs = 0
        self.__sdr = np.empty((3,))
        self.__nv = np.empty((3,))
        self.__dv = np.empty((3,))
        self.__Mm = np.empty((3,3))
        self.__Qϕδλ = np.empty((3,3))
        self.__mesh = np.zeros((3, 3))
        self.sdr_set = False
        # self.__GridAlongSD = np.empty((1,1), dtype=np.float64)
        self.__SlipGridAlongS = np.empty((1,1), dtype=np.float64)
        self.__SlipGridAlongD = np.empty((1,1), dtype=np.float64)
        
        if (strike is not None) and (dip is not None) and (rake is not None):
            print("Set Fault Segment properties:\n")
            self.SetStrikeDipRake = (strike, dip, rake)
            self.SetFaultVectors = (strike, dip, rake)
            self.SetRotationTensor = (strike, dip, rake)
            
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
    
    @property
    def RotationTensor(self):
        return self.__Qϕδλ
    
    @property
    def Mesh(self):
        """
        Return the 3D mesh for this segment

        The mesh is a 3x3 numpy array, where each row and column
        represents the x, y, and z coordinates of each node in the
        mesh. The mesh is constructed by calling the
        ``Triangulate`` method of the ``FaultMesh`` class, which
        generates a 3D mesh from the 2D grid points on the fault
        surface.

        Returns
        -------
        mesh : 3x3 numpy array
            The 3D mesh for this segment
        """
        return self.__mesh
    
    @property
    def nLsegment(self):
        return self.__nLs
    
    @property
    def nWsegment(self):
        return self.__nWs
    
    # @property
    # def GridAlongSD(self):
    #     return self.__SlipAlongSD
    
    @property
    def SlipGridAlongS(self):
        return self.__SlipGridAlongS

    @property
    def SlipGridAlongD(self):
        return self.__SlipGridAlongD

    @SlipGridAlongS.setter
    def SetSlipGridAlongS(self, GridAlongS: np.float64) -> None:
        self.__SlipGridAlongS = GridAlongS

    @SlipGridAlongD.setter
    def SetSlipGridAlongD(self, GridAlongD: np.float64) -> None:
        self.__SlipGridAlongD = GridAlongD

    # @GridAlongSD.setter
    # def SetGridAlongSD(cls, SlipFileDict: dict) -> None:
    #     cls.SlipGridAlongS = np.genfromtxt(cls.SlipFile, usecols=1)
    #     cls.SlipGridAlongD = np.genfromtxt(cls.SlipFile, usecols=0)
        
    @StrikeDipRake.setter
    def SetStrikeDipRake(self, sdr: tuple[float]) -> None:
        """
        Set Strike, Dip, Rake and convert them to radians
        
            Parameters: 
            sdr (tuple): (strike, dip, rake)
            
            Returns:
            None
        """
        strike, dip, rake = sdr
        print("Strike: {} Dip: {} Rake: {}".format(strike, dip, rake))
        self.__ϕ = strike*np.pi/180.0
        self.__δ = dip*np.pi/180.0
        self.__λ = rake*np.pi/180.0
        self.__sdr = np.array([strike, dip, rake])*np.pi/180.0
        self.sdr_set = True
    
    @NormalVector.setter
    def SetNormalVector(self, nv: list[float]) -> None:
        """
        Set fault unit-norm vector
        
            Parameters: 
            nv (list[float]): unit-norm normal vector
            
            Returns:
            None
        """
        self.__nv = nv
    
    @SlipVector.setter
    def SetSlipVector(self, dv: list[float]):
        """
        Set slip unit-norm vector
        
            Parameters: 
            dv (list[float]): unit-norm slip vector
            
            Returns:
            None
        """
        self.__dv = dv
        
    @MomentUnitTensor.setter
    def SetMomentUnitTensor(self, Mm: list[float]):
        """
        Set fault unit-norm seismic Moment
        
            Parameters: 
            Mv (list[float]): unit-norm seismic Moment
            
            Returns:
            None
        """
        self.__Mm = Mm

    @FaultVectors.setter
    def SetFaultVectors(self, sdr: tuple[float]) -> None:
        """
        Compute normal/slip vectors and unit norm Moment tensor
        
            Parameters:
            sdr (tuple): (strike, dip, rake)
            
            Returns:
            None
        """
        strike, dip, rake = sdr
        if not self.sdr_set:
            self.SetStrikeDipRake(strike, dip, rake)
        self.__ndM = compute_seismic_moment_vectors(strike=self.__ϕ,
                                                    dip=self.__δ,
                                                    rake=self.__λ)
        self.__nv, self.__dv, self.__Mm = self.__ndM
        
    @RotationTensor.setter
    def SetRotationTensor(self, sdr: tuple[float]) -> None:
        """
        Compute the fault plane rotation tensor [N,E,D] et to [E,N,Z]
        
            Parameters:            
            sdr (tuple): (strike, dip, rake)
            
            Returns:
            None
        """
        strike, dip, rake = sdr
        if not self.sdr_set:
            self.SetStrikeDipRake(strike, dip, rake)
        self.__Qϕδλ = get_rotation_tensor(strike, dip)
    
    # @Mesh.setter
    def SetMesh(self):
        """
        Set fault mesh (3D)
        
            Parameters: 
            mesh (mesh): 3D fault plane-segment mesh
            
            Returns:
            None
        """
        faultmesh = FaultMesh(xg=self.SlipGridAlongS,
                              yg=self.SlipGridAlongD,
                              zg=np.zeros_like(self.SlipGridAlongS))

        # Generate 3D flat Delaunay mesh
        faultmesh.Triangulate()
        faultmesh.RotTransMesh3d(self.__ϕ,self.__δ,self.HypoXYZ)
        
        self.__mesh = faultmesh

    @nLsegment.setter
    def SetnLsegment(self, nLs: int) -> None:
        """
        Set number of sub-faults 

        Args:
            nL (int): number of sub-faults along strike
        """
        print("Number of sub-faults along strike: {}".format(nLs))
        self.__nLs = nLs

    @nWsegment.setter
    def SetnWsegment(self, nWs: int) -> None:
        """
        Set number of sub-faults 

        Args:
            nW (int): number of sub-faults along dip
        """
        print("Number of sub-faults along dip: {}".format(nWs))
        self.__nWs = nWs