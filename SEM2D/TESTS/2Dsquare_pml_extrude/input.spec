# -*- mode: perl -*-
run_name = "Square_PML_extrude";

# duration of the run
sim_time = 2.;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=2;

snapshots {
    save_snap = true;
    snap_interval = 0.02;
};

# Description des capteurs
save_traces = true;
station_file = "capteurs.dat";
traces_format=hdf5;

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=0;

# pulse at the centre of the physical [0,500]x[0,300] domain
source {
    coords = 250. 150. 0.;
    type = impulse;
    dir = x;
    func = ricker;
    tau = .3;
    freq = 4.;
};

time_scheme {
    accel_scheme = false;
    veloc_scheme = true;
    alpha = 0.5;
    beta = 0.5;
    gamma = 1;
    courant = 0.2;
};

amortissement {
    nsolids = 0;
    atn_band = 10  0.05;
    atn_period = 0.2;
};
