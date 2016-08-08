# -*- mode: perl -*-
run_name = "Cube_SURF";

# duration of the run
sim_time = 0.008;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
mpml_atn_param = 0.002;

snapshots {
    save_snap = true;
    snap_interval = 0.001;
    deselect all;
    select box = -5  0 0  5  5 10;
    select box = -5 -5 0  0  0 10;
    select box = -5 -5 10 5  5 12;
};

# Description des capteurs
save_traces = true;
traces_format=hdf5;


# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=370;


# introduce a source
#source {
#    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
#    coords = 25. 25. 25.;
#    # the numbers before the labels are here to help convert from previous input.spec format
#    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
#    type = impulse;
#    # Direction 0.x,1.y ou 2.z (only for Impulse)
#    dir = 1. 0. 0.;
#    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
#    func = ricker;
#    tau = 0.4;
#    freq = 3.;   # source main frequency / cutoff frequency
#};

surface{
    use   = 1;              # Considérer cette surface dans le problème : 0 (Non), 1(oui)
    type  = neumann;        # type associé à la surface : neumann, planewave, fault, dirichlet
    mat_i = 0;              # indice du domaine sur lequel récupérer les propriétés matériau
    nsurf = 1;              # nombre de surfaces aux quels le type est appliqué
    index = 1;              # list de tags des surfaces concernées: ceux sont les mêmes tags donnés depuis gmsh
    C     = 0.0 0.0 0.0;    # Coordonnées du point de référence pour le calcul des différentes forces
    time  = analytic;       # fonction en temps : ricker, gauss, analytic
    freq  = 30.0;           # pour une time = ricker
    tau   = 0.0333333333;   # pour une time = ricker
    ampli = 100;            # amplitude appliquée utilisée si time est différent de analytic
    shape = paraboloid;     # formes spatiale d'une contrainte surfacique appliquée :gaussian,paraboloid,square,cylinder,uniform
    size  = 0.5;            # rayon maximun sur lequel la force s'applique si shape#uniform.
    dir   = 0.0 0.0 1.0;    # direction de la force
    var   ="t";             
    fxx   ="100.0*sign(t)";
};

surface{
    use   = 0;
    type  = dirichlet;
    mat_i = 1;
    nsurf = 9;
    index = 6 7 8 15 16 17 24 25 26;
};

surface{
    use   = 0;
    type  = planewave;
    time  = ricker;
    ampli = 500.0;
    freq  = 30.00;
    mat_i = 0;
    nsurf = 1;
    index = 1;
    dir   = 0.0  -1.0 -1.0;  # Direction de propagation de l'onde
    C     = 0. 0. 0.;
    wave  = P;               # Nature de l'onde Incidente : P (onde P), S (onde S), SH (onde SH), SV (onde SV), NO (onde non prédéfinie)
    speed = 380.;            # Vitesse de propagation de l'onde incidente : à,définir si wave=NO
    dirU  = 0.0 0.0 1.0;     # Direction du mouvement des particles : à définir si wave=NO  
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta = 0.5;            # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant = 0.2;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
};

capteurs "UU" {
    type = points;
    file = "stations.txt";
    period = 1;
};
material {
    type = constant;
};
out_variables {
    enP  = 0;   # P-wave energy (scalar field)
    enS  = 0;   # S-wave energy (scalar field)
    evol = 0;   # volumetric strain (scalar field)
    pre  = 0;   # pressure (scalar field)
    dis  = 1;   # displacement (vector field)
    vel  = 1;   #  velocity (vector field)
    acc  = 1;   # acceleration (vector field)
    edev = 0;   # deviatoric strain (tensor field)
    sdev = 0;   # deviatoric stress (tensor field)
};
