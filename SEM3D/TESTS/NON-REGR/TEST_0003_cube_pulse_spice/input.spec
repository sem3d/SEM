# -*- mode: perl -*-
run_name = "Cube_PML";

# duration of the run
sim_time = 5.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;

snapshots {
    save_snap = true;
    snap_interval = 0.02;
    select all;
};

# Description des capteurs
save_traces = true;
station_file = "capteurs.dat";

# Fichier protection reprise
prorep=false;
prorep_iter=1000;


# introduce a source
source {
    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
    coords = 0. 0. 0.;
    # the numbers before the labels are here to help convert from previous input.spec format
    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1. 0. 0.;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = spice_bench;
    tau = 0.2;
    freq = 5.;   # source main frequency / cutoff frequency
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};
