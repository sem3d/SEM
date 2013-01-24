run_name = "Toto";
accel_scheme = false;                # Acceleration scheme for Newmark
veloc_scheme = true;                 # Velocity scheme for Newmark
sim_time = 800;                    # Duration of the simulation
alpha = 0.5;                    # alpha (Newmark parameter)
beta = -0.5;                    # beta (Newmark parameter)
gamma = 1;                      # gamma (Newmark parameter)
mesh_file = "mesh4spec.";             # input mesh file
model = CUB;
anisotropy = true;                 # Model (homo, prem, 3D_Berkeley or CUB) and anisotropy (T or F)
mat_file = "elastic_n4.hetero";      # material properties file

save_snap = false;                # save snapshots
snap_interval = 100;                    # Time interval to skip snapshots

save_traces = true;                 # save traces
traces_interval = 2000;               # interval de sync des traces
station_file = "stations";          # station location file

prorep = true; # Protection reprise
prorep_iter = 1000; # Nb iter pro
verbose_level = 10; # 0, 1, 2

source {                 # introduce a source
src_coords = 29.46 65.89 6351000;    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
src_type = 2;                      # Type (1 Impulse, 2 Moment Tensor)
src_dir = 1;                      # Direction (0 x, 1 y, 2 z) (only for Impulse)
src_moment = 1e+18 1e+18 -2e+18     # M11 M22 M33 (only for Moment Tensor)
0e+18 0e+18 0e+18;      # M12 M13 M23 (only for Moment Tensor)
src_func = 3;                      # Function (1 Gaussian, 2 Ricker, 3 TF(Heaviside))
src_tau = 300;                    # tau
src_freq = 0.025;                  # source main frequency (only for Ricker)
src_band = 0.001 0.009 0.012 0.02; # source frequency band f1,f2,f3,f4 (only for TF(Heaviside))
};

neumann_cond = "Nom fichier"; # Conditions de Neumann

mpml_atn_param = 0.005; # Default 0 -> eqv PML

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
#
#amortissement {
#nsolids = 0                      # number of solids for attenuation (0 if no attenuation)
#atn_band = 1000  50;               # attenuation period band
#model_period = 1;                      # model period 
#
#};