# -*- mode: perl -*-
run_name = "Square_PML_extrude";

# duration of the run
sim_time = 2.;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=2;
ngll=5;

pml_infos {
    pml_type = PML;
};

snapshots {
    save_snap = true;
    snap_interval = 0.02;
};

# Description des capteurs -- text (.vel, legacy) AND h5 (capteurs {} below, traces_format).
# save_traces gates ALL trace output; traces_format selects the capteur format (text/hdf5).
save_traces = true;
station_file = "capteurs.dat";
traces_format=hdf5;

out_variables {
    dis = 1;
    vel = 1;
};

capteurs "REC" {
    type = points;
    file = "capteurs.dat";
    period = 1;
};

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=0;

# pulse at the centre of the physical [0,500]x[0,300] domain
source {
    coords = 250. 150.;
    type = impulse;
    dir = 1. 0.;
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
