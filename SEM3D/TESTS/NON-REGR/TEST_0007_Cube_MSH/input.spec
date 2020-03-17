# -*- mode: perl -*-
run_name = "Sm_h";

# duration of the run
sim_time =      1.500;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
mpml_atn_param=0.002;
ngll=5;

snapshots {
    save_snap = true;
    snap_interval = 0.05;
};

# Description des capteurs
save_traces = true;
traces_format = hdf5;

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=68000;

# introduce a source 0
source {
# coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
coords = 0.0  0.0  1.5;
# the numbers before the labels are here to help convert from previous input.spec format
# Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
type = impulse;
# Direction 0.x,1.y ou 2.z (only for Impulse)
dir = 0. 1. 0.;
# Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus,13.dm
func = ricker;
tau = 0.11;
freq = 10.00;
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.4;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
};

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1000;
};

material {
    type = constant;
};

out_variables {
    enP = 0;   # P-wave energy (scalar field)
    enS = 0;    # S-wave energy (scalar field)
    evol = 0;   # volumetric strain (scalar field)
    pre  = 0;   # pressure (scalar field)
    dis   = 1;   # displacement (vector field)
    vel   = 0;   #  velocity (vector field)
    acc  = 0;   # acceleration (vector field)
    edev = 0;  # deviatoric strain (tensor field)
    sdev  = 0;  # deviatoric stress (tensor field)
};
