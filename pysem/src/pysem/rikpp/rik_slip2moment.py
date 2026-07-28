# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to convert RIKsrf2 output into SEM3D input

Example:

    python3 rik_slip2moment.py @nL 140 @nW 80 @sf ./slipdistribution.dat 
    @hL 3.66 @hW 3.15 @L 7.0 @W 4.0 @fg teil @strike 45. @dip 60. @rake 108. 
    @nt 2048 @dt 0.01 @hE 631892.2 @hN 4931475.0 @hZ -860.0 @wkd ./ @tag teil
    
    
    python3 rik_slip2moment.py @nL 117 @nW 18 @sf ./slipdistribution.dat 
        @hL 32.5 @hW 17.5 @L 195.0 @W 30.0 @fg Turkey2023 @strike 276. 250.0 60.0  @dip 80.0 80.0 80.0 @rake -1.0 -1.0 -1.0 @nt 600 @dt 0.1 @hLON 37.203  @hLAT 38.024  @hZ -12960.0 @wkd ./ @tag Turkey2023-Mw7.5 @wkd /home/kltpk89/Downloads/Turkey2023/ @kf Turkey2023Mw75_kine.h5 @mf Turkey2023Mw75_moment.h5 @ct moment slip srate @snapshots 0.1 @segL 80.0 35.0 80.0 @segW 30.0 30.0 30.0 @M0 5.0484345996999996e+20


"""
# Required modules
from rik_pp_lib import *
from functools import partial
from sem_fault import SEM3Dfault, FaultSegment
from sem_fault import FaultMesh, get_segment_Lidx
from rik_fault import RIK2DFault

# General informations
__author__ = "Filippo Gatti and Elif Oral"
__copyright__ = "Copyright 2018, MSSMat UMR CNRS 8579, CentraleSupélec"
__credits__ = ["Filippo Gatti", "Elif Oral"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"




if __name__=='__main__':
    # Parse command line
    opt = start_rik()
    # Update global variables
    globals().update(opt)
    assert len(strike)==len(dip) and len(dip)==len(rake)
    
    
    fg = opj(wkd,f"{tag}")
    
    # RIK input files 
    coord_file_name = opj(wkd, 'source_coordinates.csv')
    sd_file = opj(wkd,'slipdistribution.dat')
    mr_file = opj(wkd,'MomentRate.dat')
    sr_file = opj(wkd, 'SlipRate.dat')
    hypfile = opj(wkd,'nucleationpoint.dat')

    # RIK time vector
    vtm = np.linspace(0.0, (nt-1)*dt, nt)  
    
    # Number of segments
    nSegments = len(strike)
    
    # Parse RIK 2D fault
    RIKsource = RIK2DFault(nL=int(nL), nW=int(nW))
    RIKsource.nt = int(nt)
    RIKsource.dt = dt
    RIKsource.HypoAlongSDDepthKm = {"HypoFile": hypfile,
                                    "HypoDepthKm": np.abs(hZ)/1.0e3}
    RIKsource.GridAlongSD = {"SlipFile": sd_file}
    RIKsource.MomentRateGridAlongSD = {"MomentRateFile": mr_file}
    RIKsource.SlipRateGridAlongSD = {"SlipRateFile": sr_file}
    RIKsource.HypoDepthm = 0.0
    
    # SEM3D source
    
    # Define SEM3D source
    # Initialize 
    SEMsource = SEM3Dfault(moment=M0)
    # SEM3D hypocenter (UTM) coordinates [X Y Z] = [E N Z] [m]
    SEMsource.HypoXYZ = {'ew': hE, 'ns': hN, 'ud': hZ}

    nLxSegment=list(map(partial(get_segment_Lidx, L=L, nL=nL), segL))
    nWxSegment = [nW]*nSegments
    # Initialize SEM3D fault segment left edges
    SEMsource.SegmentEdgesL = np.cumsum(np.array([0]+nLxSegment))
    # Initialize SEM3D fault segment lower edges
    SEMsource.SegmentEdgesW = np.array([0]+nWxSegment)

    # Initialize fault segments
    SEMsegments = []
    for i, (φs, δ, λ) in enumerate(zip(strike, dip, rake)):
        fs = FaultSegment(strike=φs, dip=δ, rake=λ)
        fs.SetnLsegment = nLxSegment[i]
        fs.SetnWsegment = nWxSegment[i]
        idx = np.ix_(np.linspace(SEMsource.SegmentEdgesL[i],
                                 SEMsource.SegmentEdgesL[i+1],
                                 fs.nLsegment,
                                 endpoint=False,
                                 dtype=np.int64),
                     np.linspace(SEMsource.SegmentEdgesW[0],
                                 SEMsource.SegmentEdgesW[1],
                                 fs.nWsegment,
                                 endpoint=False,
                                 dtype=np.int64)
                     )
        # Centering on the fault hypocenter
        fs.SetSlipGridAlongS = (RIKsource.SlipGridAlongS[idx]-RIKsource.HypoAlongSDDepthKm[0])*1.0e3
        fs.SetSlipGridAlongD = (RIKsource.SlipGridAlongD[idx]-RIKsource.HypoAlongSDDepthKm[1])*1.0e3
        SEMsegments.append(fs)

    print("Number of segments: {:d}".format(SEMsource.nSegments))
    
    for (i, s), d, r in zip(enumerate(strike), dip, rake):
        print(f"Segment {i+1}: strike: {s}, dip: {d}, rake: {r}\n")
            
    for i,fs in enumerate(SEMsegments):
        # Set fault-segment mesh
        nLs = fs.nLsegment
        nWs = fs.nWsegment
        nv = fs.NormalVector
        dv = fs.SlipVector
        fs.SetMesh()
        
        Moment = np.zeros((nLs, nWs, nt), dtype=np.float32)
        Slip = np.zeros((nLs, nWs, nt), dtype=np.float32)
        SlipRate = np.zeros((nLs, nWs, nt), dtype=np.float32)
        
        with h5py.File(opj(wkd,f"{tag}_kine_{i}.h5"), 'w') as kfo,\
             h5py.File(opj(wkd,f"{tag}_moment_{i}.h5"), 'w') as sfo:
        #driver='mpio', comm=MPI.COMM_WORLD)
        
            # Assign general properties
            sfo.create_dataset('vtm', data=vtm)
            kfo.attrs['dt'] = vtm[1]-vtm[0]
            kfo.attrs['Nt'] = RIKsource.nt
            kfo.attrs['Ns'] = nLs
            kfo.attrs['Nd'] = nWs
            kfo.attrs['Vnormal'] = nv
            kfo.attrs['Vslip'] = dv
            
            fs.Mesh.write_mesh2h5(sfo)
            fs.Mesh.write_mesh2h5(kfo)
            
            # Write dynamic fields        
            idx = np.ix_(np.linspace(SEMsource.SegmentEdgesL[i],
                                    SEMsource.SegmentEdgesL[i+1],
                                    nLs,
                                    endpoint=False,
                                    dtype=np.int64),
                        np.linspace(SEMsource.SegmentEdgesW[0],
                                    SEMsource.SegmentEdgesW[1],
                                    nWs,
                                    endpoint=False,
                                    dtype=np.int64),
                        np.linspace(0,
                                    nt,
                                    nt,
                                    endpoint=False,
                                    dtype=np.int64),
                        )
            
            Moment = RIKsource.MomentGridAlongSD[idx]
            Slip = RIKsource.SlipGridAlongSD[idx]
            SlipRate = RIKsource.SlipRateGridAlongSD[idx]
            
            sfo.create_dataset('moment', data=Moment, chunks=(1, 1, nt))
            
            for t in range(RIKsource.nt):
                sfo.create_dataset(f"mom_{t}", shape=(nLs*nWs,),
                                data=Moment[:, :, t].T.reshape((nLs*nWs,)),
                                chunks=(1,))
                kfo.create_dataset(f"slp_{t}", shape=(nLs*nWs,),
                                data=Slip[:, :, t].T.reshape((nLs*nWs,)),
                                chunks=(1,))
                kfo.create_dataset(f"sra_{t}", shape=(nLs*nWs,),
                                data=SlipRate[:, :, t].T.reshape((nLs*nWs,)),
                                chunks=(1,))
        # Close hdf5
        kfo.close()
        sfo.close()
        if plot:
            clr = sns.color_palette("cool",nLs*nWs)
            fig = plt.figure(figsize=(10,5))
            sns.set_style('whitegrid')
            ax  = fig.add_subplot(111)
            for i in range(0, nLs):
                for j in range(0, nWs):
                    ax.plot(vtm,Moment[i,j,:],label='Point '+str(i+1),
                            color=clr[nWs*i+j])
            # #
            ax.set_xlabel(r'$t$ [s]',fontsize=20)
            ax.set_ylabel(r'$M_0$ [Nm]',fontsize=20)
            plt.gca().spines["top"].set_alpha(0)
            plt.gca().spines["bottom"].set_alpha(.3)
            plt.gca().spines["right"].set_alpha(0)
            plt.gca().spines["left"].set_alpha(.3)
            # ax.legend()
            fig.savefig("{:>s}_moment_{:>d}_ths.png".format(fg,i), 
                        dpi=500, 
                        bbox_inches='tight')
            fig.savefig("{:>s}_moment_{:>d}_ths.eps".format(fg,i), 
                        dpi=500, 
                        bbox_inches='tight',
                    format='eps')
            plt.close()