# -*- mode: perl -*-
# 2D SOLID-FLUID interface test -- anisotropic-from-file (h5) materials.
#   solid (top, z in [0.015,0.03])  = Hooke_Aniso  <- Cstar_solid.h5 (iso-equiv vp3000 vs1700 rho2000)
#   fluid (bottom, z in [0,0.015])  = Fluid_Aniso  <- Cstar_fluid.h5 (iso-equiv vp1500 rho1000)
# Mesh: on-the-fly mesher2D, 2 layers S-over-F (mat.dat/mater.in), NOT mesh2dc.
# Source: explosion in the SOLID; a P-wave must transmit across the interface into the fluid.
run_name = "t4_solid_fluid_aniso";

sim_time = 4.0e-5;
mesh_file = "mesh4spec";               # written by mesher2D in this dir (mesh4spec.NNNN.h5)
mat_file  = "material.input";          # written by mesher2D from mater.in (S/F topology)
dim = 2;

snapshots {
    save_snap = true;
    snap_interval = 1.0e-6;
};

save_traces = true;
station_file = "stations.txt";
traces_format = hdf5;

prorep = false;
prorep_iter = 5000;
restart_iter = 0;

# explosion (isotropic moment) source, inside the SOLID top layer, off element boundaries
source {
    coords = 0.02375  0.02375 ;
    type = moment;
    moment = 1.0 1.0 0.0 ;             # Mxx Mzz Mxz (2D)
    func = ricker;
    tau = 1.5e-5;
    freq = 1.0e5;
};

time_scheme {
    accel_scheme = false;
    veloc_scheme = true;
    alpha = 0.5;
    beta  = 0.5;
    gamma = 1;
    courant = 0.3;
};

ngll = 5;

amortissement {
    nsolids = 0;
    atn_band = 10  0.05;
    atn_period = 0.2;
};

# UU_0001 solid (above source), UU_0002/0003 fluid (below interface).
# For the fluid receivers, recorded Veloc column 0 = VelPhi = -pressure.
capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1;
};

out_variables {
    dis = 1;
    vel = 1;
};
