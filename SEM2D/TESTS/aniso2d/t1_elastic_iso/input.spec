# -*- mode: perl -*-
# 2D elastic anisotropic-from-file test (isotropic-limit gate).
run_name = "t1_elastic_iso";

sim_time = 3.0e-5;
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

# explosion (isotropic moment) source, well inside the domain
source {
    coords = 0.02375  0.01875 ;
    type = moment;
    moment = 1.0 1.0 0.0 ;        # Mxx Mzz Mxz (2D)
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

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1;
};

out_variables {
    dis = 1;
    vel = 1;
};
