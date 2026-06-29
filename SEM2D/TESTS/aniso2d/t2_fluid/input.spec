# -*- mode: perl -*-
# 2D fluid anisotropic-density test (potential formulation, approach A).
# Use Cstar_fluid_iso.h5 (isotropic-limit gate) or Cstar_fluid_aniso.h5 (speed ratio).
run_name = "t2_fluid";

sim_time = 5.0e-5;
mesh_file = "sem/mesh4spec";          # generate with the 2D mesher for [0,0.05]x[0,0.03]
mat_file  = "material.input";
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

# pressure pulse in the fluid (fluidpulse = type 3; type 7/pressure is rejected in 2D)
source {
    coords = 0.025  0.015 ;
    type = fluidpulse;
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

# receiver note: for the fluid, recorded Veloc column 0 = VelPhi = -pressure.
# UU_0001 is along x and UU_0002 along z, both at distance 0.015 from the source
# (for the anisotropy speed-ratio check in test 3).
capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1;
};

out_variables {
    vel = 1;
};
