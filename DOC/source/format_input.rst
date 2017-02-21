.. -*- coding: utf-8 -*-

====================
Format de input.spec
====================

.. _input.spec:

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
sections ``source``, ``surface`` ou ``capteurs``).

Les paramètres peuvent apparaître dans un ordre quelconque au sein
d'une section (ou du corps principal), *à l'exception du paramètre*
``dim`` qui doit être **obligatoirement** déclaré avant les autres
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

  capteurs "ligne" {
    type = line;
    counti = 50;
    point0 = 0. 0. 0.;
    point1 = 0. 5000. 0.;
    periode = 1;
  };

  # champs à sortir dans les snapshots et dans les traces
  out_variables {
    enP =  0 ;
    enS =  0 ;
    evol = 1 ;
    pre =  0 ;
    dis  = 0 ;
    vel  = 0 ;
    acc  = 0 ;
    edev = 0 ;
    sdev = 0 ;
    eTotal = 0;
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
out_variables     section  pre/vel sorties    Nom de champs à sortir en output (snapshots/traces).
surface           section                     section de description des conditions appliquée à une surfaces spécifique
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
gamma             réel     --                 Paramètre pour décrire les fonctions
time_file         chaîne   --                 Fichier contenant la source ``gabor``, ``square``, ``tanh``
amplitude         réel     --                 Facteur multiplicatif appliqué à la source temporelle
Q                 réel     --                 Amplitude de la charge mobile
Y                 réel     --                 Paramètre lié au sol pour la charge mobile
X                 réel     --                 Paramètre lié au sol pour la charge mobile
L                 réel     --                 Longueur critique pour la charge mobile
v                 réel     --                 Vitesse de la charge mobile
d                 réel     --                 Distance entre les deux charges mobiles
a                 réel     --                 Distance critique entre les deux charges mobiles
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

.. _fig-source-gauss:

.. figure:: ../figures/gaussian.eps
   :scale: 60
   :align: center

- ``ricker`` :

  .. math::

     f(t) = \left(1 - 2 \left(\pi \frac{t-\tau}{T_c}\right)^2\right) \exp\left(-\left(\pi \frac{t-\tau}{T_c}\right)^2\right)

.. _fig-source-ricker:

.. figure:: ../figures/ricker.eps
   :scale: 60
   :align: center

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

.. _fig-source-heaviside:

.. figure:: ../figures/heaviside_freq.eps
   :scale: 60
   :align: center


- ``gabor`` :

  .. math::

     \sigma(t) = 2\pi f_c (t-t_0)

     f(t) = \exp(-\left(\frac{\sigma(t)}{\gamma}\right)^2) \cos(\sigma(t)+\omega) \tau

.. _fig-source-gabor:

.. figure:: ../figures/gabor_1.eps
   :scale: 60
   :align: center
.. figure:: ../figures/gabor_2.eps
   :scale: 60
   :align: center

- ``file`` : Les données sont lues dans un fichier indiqué par le paramètre ``time_file``

- ``spice_bench`` :

  .. math::

     f(t) = 1 - (1+\frac{t}{T_c})\exp(-\frac{t}{T_c})

.. _fig-source-spice_bench:

.. figure:: ../figures/spice_bench.eps
   :scale: 60
   :align: center


- ``sinus`` :

  .. math::

     f(t) = \sin(2\pi f_c (t-t_0))

.. _fig-source-sinus:

.. figure:: ../figures/sinus.eps
   :scale: 60
   :align: center

- ``square`` : Un carré *arrondi*

  .. math::

     f(t) = \frac{\exp(2.*\gamma*(x-t_0))-1.}{\exp(2.*\gamma*(x-t_0))+1}+\frac{\exp(2.*\gamma*(t_0+\tau-x))-1.}{\exp(2.*\gamma*(t_0+\tau-x))+1}

.. _fig-source-square:

.. figure:: ../figures/square.eps
   :scale: 60
   :align: center

- ``tanh``: Une tangente hyperbolique

  .. math::

     f(t) = \frac{1}{2}\tanh(\gamma*(t-t_0)+1)

.. _fig-source-tanh:

.. figure:: ../figures/tanh.eps
   :scale: 60
   :align: center


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
deselect             voir note     --                 Dé-sélection des éléments à inclure dans les snapshots
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

Va désélectionner tous les éléments, puis resélectionner tous les éléments ayant le matériau 1,
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
    counti = 50;
    point0 = 0. 0. 0.;
    point1 = 0. 5000. 0.;
    period = 1;
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

Section ``out_variables``
=========================

Chaque mot-clé est associé à un domaine spécifique requis comme sortie (snapshots / traces). Champs de sortie par défaut (lorsque tous les mots clés sont mis à 0) sont ceux de pression et de vitesse.

================  =======  =================  =================================================================
Mot-clef          Type     Valeur par défaut  Description
================  =======  =================  =================================================================
enP               bool     0                  énergie ondes P
enS               bool     0                  énergie ondes S
evol              bool     0                  déformation volumétrique
pre               bool     1                  pression
dis               bool     0                  vecteur des déplacements
vel               bool     1                  vecteur des vitesses
acc               bool     0                  vecteur des accélérations
edev              bool     0                  tenseur des déformations déviatoriques
sdev              bool     0                  tenseur des contraintes déviatoriques
eTotal            bool     0                  energies de l'ensemble des domaines (sauf PML). Energie P, Energie S, Residu PS, Energie Cinetique, Somme des Energies
================  =======  =================  =================================================================




Section ``surface``
===================

La section ``surface``  est introduite pour définir d'éventuelles conditions aux limites
imposées sur une (des) surface(s) spécifiée(s).


Les types de conditions aux limites concernés sont :

-  `condition de Neumann`

-  `condition de Dirichlet`

-  `onde plane (futur: non utilisé)`

-  `présence de PML (futur: non utilisé)`

-  `présence de faille (futur: non utilisé)`


Cette structure dans input.spec n'est utilisable que pour un maillage non structuré
importé sous le format ``.msh`` (cf les Annexes  pour plus de précision).

Déclaration d'une section `surface`
-----------------------------------

La déclaration d'un section surface dans input.spec se présente comme suit : ::

  surface {
    use   = 1;
    type  = ??;
    mat_i = 0;
    nsurf = 1;
    index = 1;
    C     = 0.0 0.0 0.0;
    time  = ricker;
    freq  = 30.0;
    tau   = 0.0333333333;
    ampli = 100;
    shape = paraboloid;
    size  = 0.5;
    dir   = 0.0 0.0 1.0;
  };

La signification des différents mots clés est consignée dans le tableau ci-dessous:

================  ========  =================  =================================================================
Mot-clef          Type      Val. par défaut    Description
================  ========  =================  =================================================================
use               entier    --                 surface activée(1) ou désactivée (0)
type              chaine    --                 condition  associée à la section
mat_i             entier    --                 domain qui fournit propriétés matériaux
nsurf             entier    --                 Nb de surfaces concerné
index             entier    --                 liste des tags de surfaces concernées
C                 réel(3)   --                 point de référence
time              chaîne    --                 type temporel appliqué
freq              réel      --                 fréquence de ricker (`time`=`ricker`)
tau               réel      --                 temps caractéristique (`time` = ricker) :math:`\tau`
ampli             réel      --                 facteur multiplicatif de la source temporelle
shape             chaîne    --                 forme spatiale de la source
size              réel      --                 rayon max de la forme spatiale
dir               réel(3)   --                 direction de la force surfacique :math:`\underline{\underline{\sigma}}.\overrightarrow{n}`
================  ========  =================  =================================================================


Valeurs des différents mots clés
--------------------------------

   #. ``type``::

        neumann, dirichlet, fault, planewave

      Les type ``planewave`` et ``fault`` ne sont pas encore disponible dans SEM

   #. ``shape``::

        gaussian, paraboloid, square, cylinder, uniform

   #. ``time``::

        ricker, gauss, analytic

Conditions de Neumann
---------------------

Conditions prédéfinies
~~~~~~~~~~~~~~~~~~~~~~

Les conditions de Neumann sont celles qui consistent à imposer des forces surfaciques
données sous la forme :math:`\underline{\underline{\sigma}}.\overrightarrow{n}`, où
:math:`\overrightarrow{n}` est la normale extérieure à la surface et `\underline{\underline{\sigma}}`
les contraintes imposées sur la surface. Il est admis dans un premier temps, que le
champs de contraintes :math:`\underline{\underline{\sigma}}` imposable sur une surface
se met sous la forme d'une fonction à varaiables séparées : une fonction `f(x,y,z)`
d'écrivant la distribution spatiale des contraintes et une fonction `g(t)` qui décrit
son allure temporelle. Dans ce cas particulier de conditions de Neumann,
les contraintes applicables sont des contraintes de direction sans cisaillement.

.. math::

   \underline{\underline{\sigma}} = \begin{bmatrix} \sigma_{11}    &    0        &     0 \\
                                                         0         & \sigma_{22} &     0 \\
                                                         0         &    0        &  \sigma_{33}
                                    \end{bmatrix}
tel que :

.. math::

   \sigma_{ii}(x,y,z,t) = k_if(x,y,z)g(t), i=1,2,3

où :math:`k_i`  sont les composantes du vecteur définissant la direction de la force
données par le mot clé ``dir``. Les seuls valeurs possibles sont 1 et 0 pour chaque
composante de ce vecteur. Dans le cas contraire, le programme redéfinit le vecteur en
attribuant la valeur 1 aux composantes non nulles

Forme spatiale (``shape``)
--------------------------

Les formes spatiales prédéfinies (les noms constituent les valeurs
que prends `shape`) sont présentées ci-dessous. Dans ces
expressions, `R` désigne la valeur attribuée à `size`, délimitant
la zone sur laquelle la forme spatiale `f(x,y,z)` est valable. Au
delà de cette zone, les contraintes sur la surfaces concernée sont
nulles. `r` désigne la distance radiale du point de coordonées
`(x,y,z)` par rapport au point de référence
:math:`C=(x_0,y_0,z_0)`.

- ``paraboloid`` :

  Cette forme est un cas non-uniforme de repartition des
  contraintes sur la surface. La répartition se limite
  essentiellement sur un disc de rayon `R`. La distribution des
  contraintes en espace présente une forme paraboloide dont
  l'expression est définie comme suit :

  .. math::

     f(x,y,z) = \sqrt{1-\frac{r^2}{R^2}}

  .. _fig-surface-para:

  .. figure:: ../figures/neu_paraboloid.png
     :scale: 60
     :align: center

- ``gaussian`` :

  IL s'agit d'une autre forme de répartition non uniforme. Cette répartition est
  une gaussienne dont l'expression est définie comme suit :

  .. math::

      f(x,y,z) = e^{-\frac{r^2}{R^2}}

  .. _fig-surface-gau:

  .. figure:: ../figures/neu_gaussian.png
     :scale: 60
     :align: center

- ``cylinder`` :

  Il s'agit d'une répartition uniforme et homogène du champ de contraintes sur un disc de rayon `R`
  centrée sur le point `C` sur la surface concernée.

  .. math::

     f(x,y,z) = 1

  .. _fig-surface-cyl:

  .. figure:: ../figures/neu_cylinder.png
     :scale: 60
     :align: center

- ``square`` :

  Il s'agit d'une répartition uniforme et homogène du champ de contraintes sur une portion carrée
  et dont la longueur de chaque côté vaut `2R` sur la surface concernée. Le carré est centré sur le point `C`.

  .. math::

     f(x,y,z) = 1

  .. _fig-source-squ:

  .. figure:: ../figures/neu_square.png
     :scale: 60
     :align: center

- ``uniform`` :

  Cette forme spatiale définie une répartition uniforme et homogène sur toute(s) la(les) surface(s)
  concernée(s).

Forme temporelle (``time``)
---------------------------

L'évolution temporelle `g(t)` prédéfinie se limite essentiellement à une forme `ricker`
et `gaussienne` déjà définie dans la section `source`.


Généralisation des conditions de Neumann
----------------------------------------

Pour appliquée une condition de Neumann ne figurant pas parmi les cas prédéfinis, SEM offre la
possibilité à l'utilisateur de définir par lui même ces contraintes impossable
:math:`\underline{\underline{\sigma}}`. Pour ce faire, il suffit de donner au mot clé
`time=analytic`. Ce qui offre la possibilité d'imposser des contraintes en cisaillement :

.. math::

   \underline{\underline{\sigma}} = \begin{bmatrix} \sigma_{11}  & \sigma_{12} & \sigma_{13} \\
                                                    \sigma_{12}  & \sigma_{22} & \sigma_{23} \\
                                                    \sigma_{13}  & \sigma_{23} & \sigma_{33}
                                    \end{bmatrix}


Ce cas nécéssite l'introduction de nouveaux mots clés dans la section `surface` pour
indiquer une imposition particulière des conditions de Neumann. Les nouveaux mots sont
consignés dans le tableau ci-dessous :

================  ========  =================  ====================================
Mot-clef          Type      Val. défaut        Description
================  ========  =================  ====================================
var               chaîne    --                 dimension en espace et temps
fxx               chaîne    --                 expression analytique
fyy               chaîne    --                 expression analytique
fzz               chaîne    --                 expression analytique
fxy               chaîne    --                 expression analytique
fxz               chaîne    --                 expression analytique
fyz               chaîne    --                 expression analytique
================  ========  =================  ====================================



#.  Le mot clé ``var`` indique les variables dont dépendent les fonctions analytiques.
    Les valeurs possibles que peut prendre ``var`` sont les suivantes : ::

     "xyzt"; "xyt"; "xzt"; "yzt"; "yt"; "zt"; "xt";
     "t"; "z"; "y"; "z"; "xyz"; "xy"; "xz"; "yz"


#. Les possibilités sont les suivantes

   -  cas pour reproduire les conditions prédéfinies lorsque la forme temporelle est
      différente de ``ricker`` et ``gauss`` ::

        surface{
           use   = 1;
           type  = neumann;
           mat_i = 0;
           nsurf = n;
           index =  ...;
           C     = 0.0 0.0 0.0;
           time  = analytic;
           shape = paraboloid;
           size  = 0.5;
           dir   = ? ? ?;
           var   ="t";
           fxx   ="g(t)";
        };

      Un example type de cette déclaration est celui ayant permis la validation. Il s'agit
      du cas-test `TEST_0007_cube_surf` que l'utilisateur peut retrouver dans les cas de
      NON-REGRESSION


   -  cas pour reproduire les conditions prédéfinies lorsque la forme temporelle est
      différente de ``ricker`` et ``gauss`` et les formes spatiales différentes de
      celles illustrées plus haut ::


        surface {
           use   = 1;
           type  = neumann;
           mat_i = 0;
           nsurf = n;
           index = ... ;
           C     = 0.0 0.0 0.0;
           time  = analytic;
           var   ="xyzt";
           fxx   ="h(x,y,z,t)";
           fyy   ="h(x,y,z,t)";
           fzz   ="h(x,y,z,t)";
        };

   -  cas permettant d'imposser en moment donc introduir des conditions en cisaillement ::


        surface {
           use   = 1;
           type  = neumann;
           mat_i = 0;
           nsurf = n;
           index = ... ;
           C     = 0.0 0.0 0.0;
           time  = analytic;
           var   ="xyzt";
           fxx   ="h(x,y,z,t)";
           fyy   ="h(x,y,z,t)";
           fzz   ="h(x,y,z,t)";
           fxy   ="h(x,y,z,t)";
           fyz   ="h(x,y,z,t)";
           fxz   ="h(x,y,z,t)";
        };


A noter :

               Dans les trois cas illustrés ci-dessus, la présence des autres mots clés reportés dans le
               premier tableau de cette section n'influence pas la condition définie. Seul les mots clés du
               deuxième tableau définissent la forme des conditions imposées.


Les fonctions analytiques sont données sous forme de chaine de caratères et construites grace à
toute une liste de fonctions analytiques élémentaires et d'opérateurs arithmétiques et aussi
des constantes prédéfinies.


Fonctions élémentaires prédéfinies
----------------------------------
==============   ================   ======================================
Mots-clé         Valeur de retour   Description
==============   ================   ======================================
cos(x)             réel             renvoie `cosinus` de x
acos(x)            réel             renvoie `arccosinus` de x
sin(x)             réel             renvoie `sinus` de x
asin(x)            réel             renvoie `arcsinus` de x
tan(x)             réel             renvoie `tangente` de x
atan(x)            réel             renvoie inverse de `tangente` de x
abs(x)             réel             renvoie  la valeur absolue de x
floor(x)           réel             renvoie la partie entière de x
exp(x)             réel             renvoie l'exponentielle  de x
log10(x)           réel             logarithme à base 10 de x
log(x)             réel             logarithme naturel de x
sqrt(x)            réel             renvoie la racine carrée de x
sinh(x)            réel             renvoie `sh` de x
cosh(x)            réel             renvoie  `ch` de x
tanh(x)            réel             renvoie `th` de x
sign(x)            réel             renvoie Heaviside de x
dirac(x)           réel             renvoie  Dirac de x
==============   ================   ======================================

Les opérateurs arithmétiques disponibles
----------------------------------------

==============   ===============================
Opérateur        Description
==============   ===============================
** ou ^          élevation à la puissance
\*               multiplication
\+               addition
\-               soustraction
/                division
(                parenthèse ouvrante
)                parenthèse fermante
==============   ===============================



Les constantes prédéfines
-------------------------

Une seule constante est prédéfinie. Il s'agit de la valeur de :math:`\pi` que l'utilisateur
peut directement utiliser en écrivant `pi`. L'utilisateur peut également définir ses propres
valeurs constantes à l'aide des mots clés consigné dans le tableau ci-dessous.

==============  ========    ===========================
Mots-clés        type       Description
==============  ========    ===========================
paramvar         entier     indicateur de présence (0/1)
npara            entier     nombre de paramètres
param            chaîne     liste des paramètres
value            réel       liste des valeurs des paramètres
==============  ========    ===========================

A noter :

        Les noms des contantes définir par l'utilisateur doivent être de même longueur que `pi`
        c'est-à-dire des nom constitués de deux caractères sans espace.


Example d'utilisation de `time = analytic`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

L'utilisateur peut utiliser l'example ci-dessous avec le cas-test `TEST_0007_cube_surf`.
Il s'agit appliquer une forme spatiale paraboloide et une évolution temporelle donnée
par `sinus`. Cet example montre également comment définir ces propres constantes ::

  surface {
      use  = 0;
      type = neumann;
      mat_i= 0;
      nsurf = 1;
      index = 1;
      C = 0.0 0.0 0.0;
      time  = analytic;
      var= "xyzt";
      fxx= "0.0";
      fyy= "0.0";
      fzz= "500.0*sqrt((1.0-(x^2+y^2+z^2)/RR^2)*sign(RR^2-(x^2+y^2
                        +z^2)))*sin(pi*t/UU)*sign(UU-t)";
      fxy= "0.0";
      fxz= "0.0";
      #
      paramvar= 1;
      npara   = 2;
      param   = "UU RR";
      value   = 0.001  0.05;
  };


Condition de Dirichlet
----------------------

Les conditions de Dirichlet traitées ici sont celle imposéer lorsqu'on fait un maillage structurée. Ces conditions consistent
à imposer une vitesse nulle et est applicable uniquement lorqu'on traite d'un problème de fluide. Une telle conditions au limite
se déclare comme suit dans ``input.spec`` : ::


  surface {
     use   = 1;
     type  = dirichlet;
     mat_i = 0;
     nsurf = n;
     index = ... ;
  };


A noter :

     L'architecture de cette section offre la possibilité d'introduite simultanément plusieurs conditions différentes.

Example :

#. dirichlet+neumann;::

     surface {
        use   = 1;
        type  = dirichlet;
        mat_i = 1;
        nsurf = m;
        index = ... ;
     };

     surface{
        use   = 1;
        type  = neumann;
        mat_i = 0;
        nsurf = n;
        index = 1 ... ;
        C     = 0.0 0.0 0.0;
        time  = analytic;
        var   ="xyzt";
        fxx   ="h(x,y,z,t)";
        fyy   ="h(x,y,z,t)";
        fzz   ="h(x,y,z,t)";
     };

#. neumann+neumann+ ... ::

     surface{
        use   = 1;
        type  = neumann;
        mat_i = 0;
        nsurf = n;
        index = 4 8 20 ... ;
        C     = 0.0 0.0 0.0;
        time  = analytic;
        var   ="xyzt";
        fxx   ="h(x,y,z,t)";
        fyy   ="h(x,y,z,t)";
        fzz   ="h(x,y,z,t)";
     };

     surface{
        use   = 1;
        type  = neumann;
        mat_i = 0;
        nsurf = m ;
        index = 10 15 ... ;
        C     = 0.0 0.0 0.0;
        time  = analytic;
        var   ="xyzt";
        fxx   ="h(x,y,z,t)";
        fyy   ="h(x,y,z,t)";
        fzz   ="h(x,y,z,t)";
        fxy   ="h(x,y,z,t)";
        fyz   ="h(x,y,z,t)";
        fxz   ="h(x,y,z,t)";
     };

     . ...; ...

 #.  Les numéros des tags fournit par le mot clé ``index`` doivent être les mêmes que
     ceux attribués aux surfaces identifiées dans ``.msh``. En présence d'un tags non
     identifié, un message d'erreur est retourné entrainant l'arrêt de tous les calculs.

Section ``pml_infos``
=====================

2 types de PML sont disponibles: splittées et convolutionnelles.

PML splittées
-------------

Mettre OPT_CPML=OFF dans cmake.

PML convolutionnelles
---------------------

Mettre OPT_CPML=ON dans cmake.

Les paramètres par défaut à ajouter dans l'input.spec sont : ::

    pml_infos {
        pml_type = CPML;
        cpml_kappa0 = 1.0;
        cpml_kappa1 = 0.0;
        cpml_rc = 0.001;
    };

``cpml_rc`` est l'impédance (taux d'absorption). ``cpml_kappa1`` permet de faire varier kappa dans la PML, en général il est mis à 0.0, mais, si le calcul est instable, jouer sur la valeur de ``cpml_kappa1`` peut améliorer la stabilité.

Note: pour que les PMLs "matchent" le domaine adjacent, il faut avoir cpml_kappa0 = 1.0 (ne pas le modifer).
