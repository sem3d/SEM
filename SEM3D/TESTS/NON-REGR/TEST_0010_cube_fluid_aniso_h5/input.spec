# -*- mode: perl -*-
run_name = "Cube_FluidAniso_H5";

# Duration: 0.15s — wave travels 225m at Vp=1500 m/s, reflections arrive >0.17s
sim_time = 0.15;
mesh_file = "mesh4spec";
mat_file = "material.input";
dim = 3;
ngll = 5;

snapshots {
    save_snap = false;
    snap_interval = 0.01;
    select all;
};

save_traces = true;
traces_format = hdf5;

capteurs "RX" {
    type = points;
    file = "stations.txt";
    period = 1;
};

# Pressure source at cube center
source {
    coords = 250. 250. 250.;
    type = fluidpulse;
    func = ricker;
    tau = 0.2;
    freq = 5.;
    amplitude = 1.e9;
};

time_scheme {
    accel_scheme = false;
    veloc_scheme = true;
    alpha = 0.5;
    beta  = 0.5;
    gamma = 1;
    courant = 0.2;
};

out_variables {
    pre = 1;   # pressure
    enP = 1;   # P-wave energy
    enK = 1;    # S-wave energy (scalar field)
    vel   = 1;   #  velocity (vector field)
    eTotal = 1; # Total energy
};
