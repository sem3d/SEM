# -*- mode: perl -*-
run_name = "Cube_PML";

# duration of the run
sim_time = 5.0000;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
mpml_atn_param = 0.002;

snapshots {
    save_snap = true;
    snap_interval = 0.01;
    deselect all;
    select box = -100 -100  100 500 500 150;
    select box = -100  100 -100 500 150 500;
    select box =  100 -100 -100 150 500 500;
};

# Description des capteurs
save_traces = true;
traces_format=hdf5;


# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=370;


# introduce a source
source {
    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
    coords = 25. 25. 25.;
    # the numbers before the labels are here to help convert from previous input.spec format
    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1. 0. 0.;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = ricker;
    tau = 0.4;
    freq = 3.;   # source main frequency / cutoff frequency
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = 0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant = 0.2;
};

ngll=5;

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 40;
};
out_variables {
    enP = 0;   # P-wave energy (scalar field)
    enS = 0;    # S-wave energy (scalar field)
    evol = 0;   # volumetric strain (scalar field)
    pre  = 0;   # pressure (scalar field)
    dis   = 1;   # displacement (vector field)
    vel   = 1;   #  velocity (vector field)
    acc  = 1;   # acceleration (vector field)
    edev = 0;  # deviatoric strain (tensor field)
    sdev  = 0;  # deviatoric stress (tensor field)
};
