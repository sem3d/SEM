# -*- mode: perl -*-
run_name = "Strain_output";

sim_time = 0.4;
mesh_file = "mesh4spec";
mat_file = "material.input";
dim = 3;
mpml_atn_param = 0.002;
ngll = 5;

snapshots {
    save_snap = false;
};

save_traces = true;
save_strain = true;
strain_period = 1;    # write one time step out of strain_period
traces_format = text;

prorep = false;

source {
    coords = 60. 70. -40.;
    type = impulse;
    dir = 1. 0. 0.;
    func = ricker;
    tau = 0.2;
    freq = 5.;
};

time_scheme {
    accel_scheme = false;
    veloc_scheme = true;
    alpha = 0.5;
    beta = 0.5;
    gamma = 1;
    courant = 0.2;
};

amortissement {
    nsolids = 0;
};

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1;
};

out_variables {
    enP = 0; enS = 0; evol = 0; pre = 0;
    dis = 1; vel = 0; acc = 0;
    edev = 0; sdev = 0;
};
