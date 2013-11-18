# -*- mode: perl -*-
run_name = "Run_3D_trial";

# duration of the run
sim_time = 5.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";

save_snap = true;
snap_interval = 0.05;

save_traces = true;
# Fichier de description des capteurs
station_file = "capteurs.dat";


#source {                 # introduce a source
## coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
#coords = 0. 0. 0.;
#type = moment;           # Type (1 Impulse, 2 Moment Tensor, fluidpulse)
#func = spice_bench;    # Function gaussian,ricker,tf_heaviside,gabor,file,spice_bench
#tau = 0.;              # tau
#moment = 0. 0. 0. 1e18 0. 0.;
#};

source {                 # introduce a source
# coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
coords = 3400. 3400. -100.;
type = pulse;
func = ricker;
tau = 0.05;
freq = 5;
dir = z;
#type = moment;           # Type (1 Impulse, 2 Moment Tensor, fluidpulse)
#func = spice_bench;    # Function gaussian,ricker,tf_heaviside,gabor,file,spice_bench
#tau = 0.0;              # tau
#moment = 1e18 1e18 1e18 0 0 0;
};


#gradient_file="gradients.dat"  # fichier gradient

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = -0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
courant = 0.05;
};

#amortissement {
#    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
#    atn_band = 10  0.05;   # attenuation period band
#    atn_period = 0.2;      # model period 
#};


#gradient {
#material = 0;
#z0 = -500.;
#z1 = 0.;
#rho_0 = 123.;
#rho_1 = 234;
#cp_0 = 1000;
#cp_1 = 2000;
#
#};
