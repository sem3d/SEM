.. -*- coding: utf-8 -*-

Description des parametres de SEM3D
===================================


Le format de fichier SEM3D a été changé pour plus de souplesse et pour
éviter des erreurs de saisie.

Le nouveau format de fichier est de la forme suivante (exemple) ::

  mot_clef = valeur;  # commentaire
  # commentaire
  section {
     mot_clef2 = "chaine";
  };
  mot_clef3 = v1 v2 v3; # un vecteur de valeurs

Les valeurs sont des entiers, des flottants, des booléens, des chaines ou des mot-clefs.

Une chaine est une suite de caractères entre guillemets (``"``), un
mot clef est une suite de caractères alphanumeriques commençant par
une lettre, et comportant des lettres, des chiffres ou le caractère
souligné (``_``).

Les sections peuvent apparaître plusieurs fois (par exemple la section ``source``).

Les paramètres peuvent apparaître dans un ordre quelconque au sein d'une section (ou du corps
principal). Un paramètre valide peut-être ignoré si il n'est pas activé par un autre paramètre :
par exemple on peut désactiver les snapshots, tout en laissant le paramètre nombre d'itération
entre snapshot.

Les mots-clef pouvant être utilisés dans le fichier (niveau 0, hors toute section) sont décrits ici :

================  =======  =================  ===========================================================
Mot-clef          type     valeur par défaut  Description
================  =======  =================  ===========================================================
run_name          chaîne   ""                 Titre de la simulation
time_scheme       section  n/a                Section de description du schéma d'intégration en temps
sim_time          réel     aucune             Durée (temps physique) de la simulation
mesh_file         chaîne   "mesh4spec"        Nom de base des fichiers maillage
mat_file          chaîne   "material.input"   Nom du fichier de description des matériaux
anisotropy        section  n/a                Description de l'anisotropie
amort             section  n/a                Description de l'amortissement
source            section  n/a                Description d'une source (peut apparaître plusieurs fois)
save_snap         bool     false              Sauvegarde des snapshots
save_interval     réel     --                 Interval (temps physique) de sauvegarde des snapshots
model
save_traces
traces_interval
station_file
prorep
prorep_iter
verbose_level
neumann_cond
================  =======  =================  ===========================================================

Description de la section ``time_scheme`` :

================  =======  =================  ===========================================================
Mot-clef          type     valeur par défaut  Description
================  =======  =================  ===========================================================
accel_scheme
veloc_scheme
alpha
beta
gamma
================  =======  =================  ===========================================================

Description de la section ``source`` :

================  =======  =================  ===========================================================
Mot-clef          type     valeur par défaut  Description
================  =======  =================  ===========================================================
coords            réel(3)  0 0 0              Position de la source
type              kw       --                 Type spatial: impulse|moment|fluidpulse
dir               kw       --                 Direction pour le type impulse ou fluidpulse (val: x|y|z)
func              kw       --                 Type temporel: gaussian|ricker|tf_heaviside|gabor|file
moment            réel(6)  --                 Moment xx yy zz xy yx xz pour le type moment
tau               réel     --                 offset de temps de démarrage de la source
freq              réel     --                 Fréquence de coupure de la fct temporelle
band              réel(4)  --                 Description des bornes pour tf_heaviside
ts                réel     --                 Pour gabor ?
time_file         chaîne   --                 Fichier contenant la source
================  =======  =================  ===========================================================


Exemple
=======

Le fichier ci-dessus correspond à celui d'un cas test ::

  # -*- mode: perl -*-
  run_name = "Run_3D_trial";
  
  # duration of the run
  sim_time = 1.0;
  mesh_file = "mesh4spec"; # input mesh file
  mat_file = "material.input";
  
  save_snap = true;
  snap_interval = 0.01;
  
  save_traces = true;
  # Fichier de description des capteurs
  station_file = "file_station";
  
  
  source {                 # introduce a source
  # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
  coords = 0. 0. 0.;
  type = impulse; # Type (1 Impulse, 2 Moment Tensor, fluidpulse)
  dir = x;                      # Direction x,y ou z (only for Impulse)
  func = ricker;          # Function gaussian,ricker,tf_heaviside,gabor,file
  tau = 0.2;              # tau
  freq = 5.;           # source main frequency (only for Ricker)
  };
  
  
  #gradient_file="gradients.dat"  # fichier gradient
  
  time_scheme {
      accel_scheme = false;  # Acceleration scheme for Newmark
      veloc_scheme = true;   # Velocity scheme for Newmark
      alpha = 0.5;           # alpha (Newmark parameter)
      beta = -0.5;           # beta (Newmark parameter)
      gamma = 1;             # gamma (Newmark parameter)
  };
