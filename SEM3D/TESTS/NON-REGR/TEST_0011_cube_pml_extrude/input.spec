# -*- mode: perl -*-
run_name = "Cube_PML_extrude";

# duration of the run
sim_time = 2.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
mpml_atn_param = 0.002;

snapshots {
    save_snap = true;
    snap_interval = 0.02;
};

# Description des capteurs
save_traces = true;
traces_format=hdf5;

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=0;

# a pulse at the centre of the physical (non-PML) cube [0,300]^3
source {
    coords = 150. 150. 150.;
    type = impulse;
    dir = 1. 0. 0.;
    func = ricker;
    tau = 0.3;
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

ngll=5;

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 20;
};
out_variables {
    dis = 1;
    vel = 1;
};
