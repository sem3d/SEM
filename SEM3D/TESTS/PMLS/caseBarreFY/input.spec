# -*- mode: perl -*-
run_name = "Test_name";

# duration of the run
sim_time = 5.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
fmax=1.0;

snapshots {
    save_snap = true;
    snap_interval = 0.1;
    select all;
};

# Description des capteurs
save_traces = true;
traces_format=text;

capteurs "A0" {
    type = single;
    point0 = 0. 2000. 0.;
    period = 1;
};
capteurs "A1" {
    type = single;
    point0 = 0. -2000. 0.;
    period = 1;
};
capteurs "B0" {
    type = single;
    point0 = 0. 4000. 0.;
    period = 1;
};
capteurs "B1" {
    type = single;
    point0 = 0. -4000. 0.;
    period = 1;
};
capteurs "C0" {
    type = single;
    point0 = 0. 6000. 0.;
    period = 1;
};
capteurs "C1" {
    type = single;
    point0 = 0. -6000. 0.;
    period = 1;
};
capteurs "D0" {
    type = single;
    point0 = 0. 8000. 0.;
    period = 1;
};
capteurs "D1" {
    type = single;
    point0 = 0. -8000. 0.;
    period = 1;
};
capteurs "E0" {
    type = single;
    point0 = 0. 10000. 0.;
    period = 1;
};
capteurs "E1" {
    type = single;
    point0 = 0. -10000. 0.;
    period = 1;
};

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=0;

# introduce a source
source {
    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
    coords = 0. 0. 0.;
    # the numbers before the labels are here to help convert from previous input.spec format
    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
    type = fluidpulse;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = ricker;
    tau = 0.4;
    freq = 3.;
    amplitude = 1000000000.;
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta =  0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};

pml_infos {
pml_type = CPML;
cpml_kappa0 = 1.0;
cpml_kappa1 = 0.0;
cpml_rc = 0.000001;
};

