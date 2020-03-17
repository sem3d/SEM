# -*- mode: perl -*-
run_name = "%(name)s";

# duration of the run
sim_time = %(sim_time)f;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
ngll=5;
snapshots {
    save_snap = %(save_snap)s;
    snap_interval = %(snap_interval)s;
    deselect all;
    select box = %(xmin)s %(ymin)s  %(z0)s %(xmax)s %(ymax)s %(z1)s;
    select box = %(xmin)s %(y0)s  %(zmin)s %(xmax)s %(y1)s %(zmax)s;
    select box = %(x0)s %(ymin)s  %(zmin)s %(x1)s %(ymax)s %(zmax)s;
};

# Description des capteurs
save_traces = true;
station_file = "capteurs.dat";
traces_format=hdf5;


# Fichier protection reprise
prorep=%(prorep)s;
prorep_iter=1000;
restart_iter=%(restart_iter)s;


# introduce a source
source {
    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
    coords = %(xs)s %(ys)s %(zs)s;
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
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};

amortissement {
    nsolids = %(nsolids)d;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};
