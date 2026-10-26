# -*- mode: perl -*-
run_name = "Source_fluid_solid_interface";

sim_time = 0.6;
mesh_file = "mesh4spec";
mat_file = "material.input";
dim = 3;
mpml_atn_param = 0.002;
ngll = 5;

snapshots {
    save_snap = false;
};

save_traces = true;
traces_format = text;

prorep = false;

source {
    coords = 213. 237. -200.;
    type = impulse;
    dir = 1. 0. 0.;
    func = ricker;
    tau = 0.3;
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
