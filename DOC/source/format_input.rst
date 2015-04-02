.. -*- coding: utf-8 -*-

====================
Format de input.spec
====================


Syntaxe du fichier
==================

Le format de fichier :file:`input.spec` a été changé depuis la version 2012, pour plus de
souplesse et pour éviter des erreurs de saisie.

Le nouveau format de fichier est de la forme suivante (exemple) ::

  mot_clef = valeur;  # commentaire
  # commentaire
  section {
     mot_clef2 = "chaine";
  };
  mot_clef3 = v1 v2 v3; # un vecteur de valeurs

Les valeurs sont des entiers, des flottants, des booléens, des chaînes
ou des mot-clefs.

Une chaîne est une suite de caractères entre guillemets (``"``), un
mot clef est une suite de caractères alphanumériques commençant par
une lettre, et comportant des lettres, des chiffres ou le caractère
souligné (``_``).

Certaines sections peuvent apparaître plusieurs fois (par exemple les
sections ``source`` ou ``capteurs``).

Les paramètres peuvent apparaître dans un ordre quelconque au sein
d'une section (ou du corps principal), *à l'exception du paramètre*
``dim`` qui doit être **obligatoirement** declaré avant les autres
pour valider les tailles des vecteurs.

Un paramètre valide peut-être ignoré si il n'est pas activé par un
autre paramètre : par exemple on peut désactiver les snapshots, tout en
laissant le paramètre nombre d'itérations entre snapshot.

Exemple
-------

Le fichier suivant correspond à celui d'un cas test : ::

  # -*- mode: perl -*-
  run_name = "Run_3D_trial";
  
  # duration of the run
  sim_time = 1.0;
  mesh_file = "mesh4spec"; # input mesh file
  mat_file = "material.input";
  dim = 3;
  
  snapshots {
    save_snap = true;
    snap_interval = 0.01;
    deselect all;
    select material = 1;
    select box = -10 -10 -10 10 10 10;
  };
  save_traces = true;

  capteurs "AB" {
    type = points;
    # Fichier de description des capteurs
    file = "file_station";
    period = 2; # Une iteration sur 2
  };
  
  
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
    beta  = 0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
  };


Section globale
===============

Les mots-clef pouvant être utilisés dans le fichier (niveau 0, hors toute section) sont décrits ici :

.. tabularcolumns:: |p{2.8cm}|p{1.5cm}|p{2.3cm}|p{8cm}|

================  =======  =================  ================================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ================================================================
amortissement     section  n/a                Description de l'amortissement
fmax              réel     1Hz                Fréquence max attendue du signal (utilisé pour vérifications)
ngll              entier   5                  (futur) nombre de points de gauss par maille
dim               entier   obligatoire        Spécifie si le calcul est 2D ou 3D.
mat_file          chaîne   "material.input"   Nom du fichier de description des matériaux
mesh_file         chaîne   "mesh4spec"        Nom de base des fichiers maillage
mpml_atn_param    réel     0.0                Coefficient d'amortissement MPML (et activation MPML si non nul)
prorep            bool     false              Reprise d'un calcul précédent
prorep_iter       entier   n/a                Interval entre 2 protections (ou 0 pour désactiver)
restart_iter      entier   n/a                Numéro de la protection pour reprendre le calcul
run_name          chaîne   ""                 Titre de la simulation
snapshots         section  n/a                Description des paramètres de sauvegarde des snapshots
save_traces       bool     false              Activation des capteurs
sim_time          réel     aucune             Durée (temps physique) de la simulation
source            section  n/a                Description d'une source (peut apparaître plusieurs fois)
traces_interval   entier                      Interval de sortie des capteurs en nombre d'itérations
traces_format     kw       text               Format des sorties capteurs ``text`` ou ``hdf5``
time_scheme       section  n/a                Section de description du schéma d'intégration en temps
pml_info          section                     Pour l'instant 2D seul. Description des PMLs
anisotropy        bool                        (futur: non utilisé)
gradient          section                     (futur: non utilisé)
model             section                     (futur: non utilisé)
neumann           section                     (futur: non utilisé)
verbose_level     entier
capteurs          section                     Définition d'un ensemble de capteurs
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
================  ========  =================  ===========================================================

Section ``amortissement``
=========================

================  =======  =================  ===========================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ===========================================================
nsolids           entier   0                  Nombre de mécanismes. 0 signifie désactivation.
atn_band          réel(2)  n/a                Période max et min à atténuer
atn_period        réel     n/a                Période de définition de Qp et Qs
================  =======  =================  ===========================================================

Paramétrage de l'atténuation
----------------------------

Le mécanisme d'atténuation est décrit en deux endroits :

- Le fichier de description des matériaux contient les paramètres :math:`Q_\kappa` et :math:`Q_\mu` du
  milieu.

- Le fichier ``input.spec`` contient la section ``amortissement`` décrite ci-dessus.

L'atténuation est modélisée par N mécanismes (``nsolids``) sur une bande
de fréquences décrite par ``atn_band``. Les N mécanismes sont centrés sur
N fréquences choisies dans la bande spécifiée. 

Le paramètre ``atn_period`` spécifie la période pour laquelle les
valeurs de :math:`V_p` et :math:`V_s` sont exactement celles
spécifiées dans le fichier matériau.

Le code n'applique pas d'atténuation si ``nsolids=0``.


Section ``time_scheme``
=======================

.. tabularcolumns:: |p{3cm}|p{1.5cm}|p{2cm}|p{8cm}|

================  =======  =================  ===========================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  ===========================================================
accel_scheme      bool                        Schéma en temps
veloc_scheme      bool                        Schéma en vitesse
alpha             réel                        Paramètre :math:`\alpha` d'intégration de Newmark
beta              réel                        Paramètre :math:`\beta` d'intégration de Newmark
gamma             réel                        Paramètre :math:`\gamma` d'intégration de Newmark
courant           réel     0.2                Nombre de Courant. Le calcul du pas de temps en dépend.
================  =======  =================  ===========================================================

Section ``source``
==================

================  =======  =================  =================================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  =================================================================
coords            réel(3)  0 0 0              Position de la source
type              kw       --                 Type spatial: impulse|moment|fluidpulse
dir               kw       --                 Direction pour le type impulse ou fluidpulse (val: x|y|z)
func              kw       --                 Type temporel (voir :ref:`defsources` ci-dessous)
moment            réel(6)  --                 Moment xx yy zz xy xz yz pour le type moment
tau               réel     --                 Un temps caractéristique :math:`\tau`
freq              réel     --                 Une fréquence :math:`f_c`
band              réel(4)  --                 Description des bornes :math:`f_1,f_2,f_3,f_4` pour tf_heaviside
ts                réel     --                 Un offset de temps :math:`t_0`
gamma             réel     --
time_file         chaîne   --                 Fichier contenant la source
amplitude         réel     --                 Facteur multiplicatif appliqué à la source temporelle
Q                 réel     --                 
Y                 réel     --                 
X                 réel     --                  
L                 réel     --                 
v                 réel     --                 
d                 réel     --                 
a                 réel     --                 
================  =======  =================  =================================================================

Note:
  Depuis la version 2014.09, la dimension des vecteurs et matrices ci-dessus, dépend de la dimension
  du problème (paramètre dim=2 ou dim=3). En 2D les paramètres ``coords`` et ``moment`` sont respectivement
  de dimension 2 et 4.

.. _defsources:

Paramètres des sources
----------------------

Les formes d'ondes temporelles des sources sont décrites ci-dessous. Les
paramètres sont décrits dans la section ``source``. Certains sont calculés :

  - :math:`f_c` : paramètre ``freq``
  
  - :math:`T_c = \frac{1}{f_c}`
  
  - :math:`\tau` : paramètre ``tau``
  
  - :math:`t_0` : paramètre ``ts``
  
  - :math:`f_1,f_2,f_3,f_4` : décrits par le paramètre (4 composantes) ``band``
  
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

     f(t) = \sin(2\pi f_c (t-t_0))


- ``square`` : Un carré *arrondi*
 
  .. math::

     f(t) = \frac{\exp(2.*\gamma*(x-t_0))-1.}{\exp(2.*\gamma*(x-t_0))+1}+\frac{\exp(2.*\gamma*(t_0+\tau-x))-1.}{\exp(2.*\gamma*(t_0+\tau-x))+1}

- ``tanh``: Une tangente hyperbolique

  .. math::

     f(t) = \frac{1}{2}\tanh(\gamma*(t-t_0)+1)

- ``dm``: M function

  .. math::

     f(t) = \frac{Q*Y}{2}*(X^{\frac{v*(t-t_0)-a}{d^2}}+X^{\frac{v*(t-t_0)-a-L}{d^2}}) 


Section ``snapshots``
=====================

.. tabularcolumns:: |p{3cm}|p{1.5cm}|p{1.5cm}|p{8cm}|

===================  ============  =================  ============================================================
Mot-clef             Type          Valeur par défaut  Description
===================  ============  =================  ============================================================
save_snap            bool          false              Sauvegarde des snapshots
save_interval        réel          --                 Interval (temps physique) de sauvegarde des snapshots
select               voir note     --                 Sélection des éléments à inclure dans les snapshots
deselect             voir note     --                 Désélection des éléments à inclure dans les snapshots
group_outputs        entier        32                 Écriture d'un fichier sortie par *group_outputs* processeurs
output_total_energy  bool                             2D uniquement, calcul de l'énergie totale
===================  ============  =================  ============================================================

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



Section ``capteurs``
====================

Le mot-clef ``capteurs`` doit être suivi d'une chaîne avant le début de la section, comme dans les exemples
suivants ::

  capteurs "VERT_PT_A" {
    type = points;
    file = "cpt_vert_a.txt";
  };

  capteurs "ligne" {
    type = line;
    count = 50;
    point0 = 0. 0. 0.;
    point1 = 0. 5000. 0.;
    periode = 1;
  };

Les mots-clefs admissibles dans une section ``capteurs`` sont :

===================  ============  =================  ========================================================
Mot-clef             Type          Valeur par défaut  Description
===================  ============  =================  ========================================================
type                 kw            --                 Type de description (point,line,plane,single)
file                 fichier       --                 Chemin vers un fichier
counti               entier        --                 Nombre de points (:math:`N_i`)
countj               entier        --                 Nombre de points (:math:`N_j`)
period               entier        1                  Période de sortie du (groupe de) capteur(s)
point0               coordonnées   --                 Point 0 (:math:`P_0`)
point1               coordonnées   --                 Point 1 (:math:`P_1`)
point2               coordonnées   --                 Point 2 (:math:`P_2`)
===================  ============  =================  ========================================================

Description des type de capteurs :

- ``points`` : Une liste de points, définis dans un fichier spécifié par le mot-clef ``file``

- ``single`` : Le plus simple, défini par le mot clef ``point0``, son nom sera le nom de la section.

- ``line`` : Définit :math:`N_i` capteurs sur le segment :math:`[P_0, P_1]`. Un numéro leur est attribué
  leur nom est préfixé du nom de la section

- ``plane`` : Définit :math:`N_i \times N_j` capteurs sur le
  parallélépipède définit par les deux vecteurs :math:`\overrightarrow{P_0{}P_1}`
  et :math:`\overrightarrow{P_0{}P_2}`.  Un numéro leur est attribué. Leur nom est
  préfixé du nom de la section.

  Pour ``i`` variant de 0 à :math:`N_i-1` et ``j`` variant de 0 à
  :math:`N_j-1`, alors les coordonnées des points sont : :math:`P_{ij}
  = P_0 + \frac{i}{N_i-1} \overrightarrow{P_0 P_1} + \frac{j}{N_j-1}
  \overrightarrow{P_0 P_2}`


