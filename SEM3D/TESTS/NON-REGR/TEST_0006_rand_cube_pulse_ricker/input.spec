# -*- mode: perl -*-
run_name = "Cube_RAND";

# duration of the run
sim_time =      5.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
mpml_atn_param=0.002;

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
prorep_iter=4000;
restart_iter=00452000;


# introduce a source 0
source {
# coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
coords =        0.00.        0.00.        0.00.;
# the numbers before the labels are here to help convert from previous input.spec format
# Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
type = impulse;
# Direction 0.x,1.y ou 2.z (only for Impulse)
dir = 0. 1. 0.;
# Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus,13.dm
func = ricker;
tau = 0.0333333333;
freq =       30.00;
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.12;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period
};
capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 40;
};
material {
    type = random;
    random_library_path = "/home/carvalhol/Projects/SEM/build_RF";
};

out_variables {
    enP = 1;   # P-wave energy (scalar field)
    enS = 1;    # S-wave energy (scalar field)
    evol = 0;   # volumetric strain (scalar field)
    pre  = 0;   # pressure (scalar field)
    dis   = 1;   # displacement (vector field)
    vel   = 0;   #  velocity (vector field)
    acc  = 0;   # acceleration (vector field)
    edev = 0;  # deviatoric strain (tensor field)
    sdev  = 0;  # deviatoric stress (tensor field)
};
