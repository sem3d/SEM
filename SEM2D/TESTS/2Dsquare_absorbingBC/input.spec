# -*- mode: perl -*-
run_name = "Square_Absorbing_BC";

# duration of the run
sim_time = 0.01;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";

snapshots {
    save_snap = true;
    snap_interval = 0.00005;
    deselect all;
    select box = -100 -100  100 500 500 150;
    select box = -100  100 -100 500 150 500;
    select box =  100 -100 -100 150 500 500;
};

# Description des capteurs
save_traces = true;
station_file = "capteurs.dat";
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
    dir = x;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = ricker;
    tau = 0.4;
    freq = 3.;   # source main frequency / cutoff frequency
};

time_scheme {
    type_time_integration = 2;  # Type of time integration (1 for Newmark, 2 for RK4)
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};

type_elements {
    dg_type = 0;           # Type of DG (0 for DG-strong, 1 for DG-weak, 2 for CG)
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};

