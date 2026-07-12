# -*- mode: perl -*-
run_name = "Square_Absorbing_BC";

# duration of the run
sim_time = 2.;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim = 2;
ngll = 5;   # NGLL commun a tout le domaine (lu depuis input.spec, comme en 3D)

snapshots {
    save_snap = true;
    snap_interval = 0.02;
    deselect all;
    select box = -100 -100  100 500 500 150;
    select box = -100  100 -100 500 150 500;
    select box =  100 -100 -100 150 500 500;
};

# Description des capteurs -- h5 (capteurs {} below, traces_format). Legacy .vel writer is
# disabled in the solver; save_traces still gates all trace output.
# STATUS (2026-07-12): input.spec parse bugs fixed (added dim=2; 2D source coords/dir). The
# example still does NOT run end-to-end: (a) the shipped mesh4spec.0000 is an obsolete ASCII
# format -- the current solver needs an HDF5 mesh (regenerate via mesher2D with current-format
# mat.dat/mater.in/mesh.input); (b) with a regenerated mesh, the DG path (dg_strong/godunov +
# absorbing) crashes at runtime -- a separate DG-solver bug, out of scope here.
save_traces = true;
station_file = "capteurs.dat";
traces_format=hdf5;

out_variables {
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
restart_iter=370;


# introduce a source
source {
    # 2D coordinates (x z); domain is [-100,500] x [-100,500]
    coords = 0. 0.;
    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
    type = impulse;
    # Direction as a 2D vector (x z), like the other SEM2D decks
    dir = 1. 0.;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = ricker;
    tau = .5;
    freq = 3.;   # source main frequency / cutoff frequency
};

time_scheme {
    type_time_integration = RK4;  # Type of time integration (1 for Newmark, 2 for RK4)
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.4;
};

type_elements {
    dg_type = dg_strong;           # Type of DG (0.continuous, 1.dg_strong, 2.dg_weak)
    flux_type = godunov;     	   # Type of Flux (0.none, 1.centered, 2.godunov, 3.laurent)
    bc_type = absorbing;	   # Type of Boundary Conditions (0.free, 1.absorbing)
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};

