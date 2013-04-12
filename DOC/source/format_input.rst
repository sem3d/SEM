.. -*- coding: utf-8 -*-

===================================
Description des parametres de SEM3D
===================================

Format
======

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

================  =======  =================  ================================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ================================================================
amortissement     section  n/a                Description de l'amortissement
mat_file          chaîne   "material.input"   Nom du fichier de description des matériaux
mesh_file         chaîne   "mesh4spec"        Nom de base des fichiers maillage
mpml_atn_param    réel     0.0                Coéfficient d'amortissement MPML (et activation MPML si non nul)
prorep            bool     false              Reprise d'un calcul précédent
prorep_iter       entier   n/a                Numéro de la protection pour reprendre le calcul
run_name          chaîne   ""                 Titre de la simulation
snapshots         section  n/a                Description des paramètres de sauvegarde des snapshots
save_traces       bool     false              Activation des capteurs
traces_format     kw       text               Format des sorties capteurs ``text`` ou ``hdf5``
sim_time          réel     aucune             Durée (temps physique) de la simulation
source            section  n/a                Description d'une source (peut apparaître plusieurs fois)
station_file      chaîne   "capteurs.dat"     Fichier de description des capteurs
time_scheme       section  n/a                Section de description du schéma d'intégration en temps
verbose_level     entier
================  =======  =================  ================================================================

Les paramètres suivants sont reconnus mais non utilisés dans cette version :

================  ========  =================  ===========================================================
Mot-clef          Type      Valeur par défaut  Description
================  ========  =================  ===========================================================
anisotropy        bool      n/a                Description de l'anisotropie
gradient          section   n/a                Description des gradients
model             kw        --                 CUB|homo|prem|3D_berkeley
neumann           bool                         .
traces_interval   entier                       .
traces_format     kw
================  ========  =================  ===========================================================

Description de la section ``amortissement`` :

================  =======  =================  ===========================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ===========================================================
nsolids           entier   0                  Nombre de mécanisme. 0 signifie désactivation.
atn_band          réel(2)  n/a                Période max et min à atténuer
atn_period        réel     n/a                Période centrale (un mécanisme aura cette valeur centrale)
================  =======  =================  ===========================================================

Description de la section ``time_scheme`` :

================  =======  =================  ===========================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ===========================================================
accel_scheme      bool                        Schéma en temps
veloc_scheme      bool                        Schéma en vitesse
alpha             réel                        Paramètre :math:`\alpha` d'intégration de Newmark
beta              réel                        Paramètre :math:`\beta` d'intégration de Newmark
gamma             réel                        Paramètre :math:`\gamma` d'intégration de Newmark
courant           réel     0.2                Nombre de courant. Le calcul du pas de temps en dépend.
================  =======  =================  ===========================================================

Description de la section ``source`` :

================  =======  =================  =================================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  =================================================================
coords            réel(3)  0 0 0              Position de la source
type              kw       --                 Type spatial: impulse|moment|fluidpulse
dir               kw       --                 Direction pour le type impulse ou fluidpulse (val: x|y|z)
func              kw       --                 Type temporel (voir `Les sources`_ ci-dessous)
moment            réel(6)  --                 Moment xx yy zz xy yx xz pour le type moment
tau               réel     --                 Un temps caractéristique :math:`\tau`
freq              réel     --                 Une fréquence :math:`f_c`
band              réel(4)  --                 Description des bornes :math:`f_1,f_2,f_3,f_4` pour tf_heaviside
ts                réel     --                 Un offset de temps :math:`t_0`
gamma             réel     --
time_file         chaîne   --                 Fichier contenant la source
amplitude         réel     --
================  =======  =================  =================================================================

Description de la section ``snapshots`` :

===============  ============  =================  ============================================================
Mot-clef         Type          Valeur par défaut  Description
===============  ============  =================  ============================================================
save_snap        bool          false              Sauvegarde des snapshots
save_interval    réel          --                 Interval (temps physique) de sauvegarde des snapshots
select           voir note     --                 Sélection des éléments à inclure dans les snapshots
deselect         voir note     --                 Désélection des éléments à inclure dans les snapshots
===============  ============  =================  ============================================================

Note:
  Par défaut, les snapshots incluent toutes les mailles. Le format de la commande select/deselect
  est décrit ci-dessous.

On peut choisir de sélectionner ou déselectionner des mailles pour les inclure ou les exclure des sorties.

Il y a pour l'instant deux critères de sélection : le numéro du matériau ou la localisation absolue.

Les commandes de sélection/déselection sont appliquées dans l'ordre du fichier ``input.spec``.

La syntaxe de la commande est : ::

  [de]select (all|material = NN|box = x0 y0 z0 x1 y1 z1) ;

Ainsi : ::

  deselect all;
  select material = 1;
  selec box = -500 -10 -10 500 10 10;

Va déselectionner tous les éléments, puis resélectionner tous les éléments ayant le matériau 1,
ainsi que tous les éléments dont le centre se situe dans la boite spécifiée.

Autre exemple : ::

  select all;  # Inutile car par défaut
  deselect material  = 5;
  deselect material  = 6;
  deselect material  = 7;

Cette description va simplement exclure les matériaux 5, 6 et 7 des sorties.

Exemple
=======

Le fichier suivant correspond à celui d'un cas test : ::

  # -*- mode: perl -*-
  run_name = "Run_3D_trial";
  
  # duration of the run
  sim_time = 1.0;
  mesh_file = "mesh4spec"; # input mesh file
  mat_file = "material.input";
  
  snapshots {
    save_snap = true;
    snap_interval = 0.01;
    deselect all;
    select material = 1;
    select box = -10 -10 -10 10 10 10;
  };
  save_traces = true;
  # Fichier de description des capteurs
  station_file = "file_station";
  
  
  source {                 # introduce a source
    # coordinates of the sources
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


Les sources
===========

Les formes d'ondes temporelles des sources sont décrites ci-dessous. Les
paramètres sont décrits dans la section ``source``. Certains sont calculés :

  - :math:`f_c` : paramètre ``freq``
  
  - :math:`T_c = \frac{1}{f_c}`
  
  - :math:`\tau` : paramètre ``tau``
  
  - :math:`t_0` : paramètre ``ts``
  
  - :math:`f_1,f_2,f_3,f_4` : décrit par le paramètre (4 composantes) ``band``
  
  - :math:`\gamma` : paramètre ``gamma``


Les fonctions temporelles sont:

- ``gaussian`` :  

  .. math::

     f(t) = -2 (t-t_0) \exp\left(-\frac{(t-t_0)^2}{\tau^2}\right)

- ``ricker`` :

  .. math::

     f(t) = \left(1 - 2 \left(\pi \frac{t-\tau}{T_c}\right)^2\right) \exp\left(-\left(\pi \frac{t-\tau}{T_c}\right)^2\right)

- ``tf_heaviside`` :

  .. math::
     :nowrap:

     \begin{eqnarray}
     f(t) & = & \mathcal{TF}^{-1}(\phi(\omega)) \\
     \phi(\omega) & = & \exp(-i\omega\tau).\chi_{f_1,f_2,f_3,f_4}(\frac{\omega}{2\pi}) \\
     \chi(f) & = & 1 \text{ if } f_2 < f < f_3 \\
             &   & 0 \text{ if } f  < f_1 \text{ or } f > f_4 \\
             &   & \frac{1}{2}\left(1+\cos\left(\pi\frac{f-f_3}{f_4-f_3}\right)\right) \text{ if } f_3 < f < f_4 \\
             &   & \frac{1}{2}\left(1+\cos\left(\pi\frac{f-f_2}{f_2-f_1}\right)\right) \text{ if } f_1 < f < f_2
     \end{eqnarray}

- ``gabor`` :

  .. math::

     \sigma(t) = 2\pi f_c (t-t_0)

     f(t) = \exp(-\left(\frac{\sigma(t)}{\gamma}\right)^2) \cos(\sigma(t)+\omega) \tau

- ``file`` : Les données sont lues dans un fichier indiqué par le paramètre ``time_file``

- ``spice_bench`` :

  .. math::

     f(t) = 1 - (1+\frac{t}{T_c})\exp(-\frac{t}{T_c})

- ``sinus`` :

  .. math::

     f(t) = sin(2\pi f_c (t-t_0))

