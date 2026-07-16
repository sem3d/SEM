.. -*- coding: utf-8 -*-

============================
Format des fichiers matériau
============================

.. _material.input: 

Format de mat.dat
=================

Ce fichier est seulement nécessaire dans le cas où le maillage est réalisé par le mailleur automatique.
Le fichier ``mat.dat`` doit contenir (les commentaires, après le *#*
sont facultatifs) ::

  -100.    # xmin
  500.     # xmax
  50.      # xstep
  -100.    # ymin
  500.     # ymax
  50.      # ystep
  500.     # zmax
  2        # nb. of layers : i-th layer associated to i-th material (defined in mater.in)
  600. 12  # upper layer: thickness and nb of steps       ...
             ... from zmax, 12 elements are created backwards over 600. m => zmin = -100
  100. 1   # lower layer: thickness and nb of steps       ...
             ... from upper layer, 1 element is created backwards over 100. m => zmin = -200
  1        # PMLs on X / Y ? 0: no, 1: yes => PML on X and Y (at left and rigth for each direction)
  1 1      # PMLs on Z     ? on top ? at the bottom ? (0: no, 1: yes)
  1        # 8 or 27 control points for elements (1 or 2)

Note : après la ligne qui définit le nombre de couches, il doit y avoir autant de lignes que de couches.

- Choix de 8 noeuds par maille : 1 (Les mailles quadratiques à 27
  noeuds sont en développement)

Ce format à été étendu pour rendre plus flexible la définition des PMLs :

- Le premier indicateur (PML on X/Y), indique maintenant le nombre de couches PML : donc une valeur
  de 0 ou 1 reste compatible (présence ou non d'une couche PML), une valeur de 3 indiquera qu'on souhaite
  3 couches de PML.

- Le sens des deux indicateurs suivants (PML on top PML at the bottom)
  est modifié si on indique un 2 pour le premier indicateur. Dans ce
  cas il faut fournir non pas 1 mais 6 valeurs à la suite du 2, une pour chaque direction
  de PML dans l'ordre Up  (Z+), Down (Z-), Nord (Y+) Sud (Y-) Est (X+) West (X-)

Voici un second exemple de fichier mat.dat qui spécifie 3 couches de PML dans les directions X+, X- et Z+ ::

  -100.    # xmin
  500.     # xmax
  50.      # xstep
  -100.    # ymin
  500.     # ymax
  50.      # ystep
  500.     # zmax
  2        # nb. of layers : i-th layer associated to i-th material (defined in mater.in)
  600. 12  # upper layer: thickness and nb of steps       ...
             ... from zmax, 12 elements are created backwards over 600. m => zmin = -100
  100. 1   # lower layer: thickness and nb of steps       ...
             ... from upper layer, 1 element is created backwards over 100. m => zmin = -200
  3        # Number of extra PML layers (0=no pmls)
  2 1 0 0 0 1 1   # PMLs on U/D/N/S/E/W
  1        # 8 or 27 control points for elements (1 or 2)
  

Format de mater.in
==================

**Format commun à SEM2D et SEM3D.**

Le fichier ``mater.in`` décrit combien de matériaux sont utilisés dans le modèle ::

  1
  S  6300.00  2500.00   2800. 630. 250.

 `1` est le nombre de matériaux dans le modèle.

La deuxième ligne décrit le type de matériau (``S`` matériau solide et
``F`` matériau fluide). Pour chaque matériau, on déclare
successivement, la vitesse de propagation de l'onde de pression,
vitesse de l'onde de cisaillement, la densité du matériau,
et les paramètres :math:`Q_\kappa` et :math:`Q_\mu`
pour l'atténuation des ondes P et S.

Le nombre de points de Gauss (NGLL) n'est **pas** indiqué ici : il est commun à
tout le domaine et provient de ``input.spec`` (paramètre ``ngll=``), en 2D comme
en 3D. Le pas de temps n'y figure pas non plus : il est calculé à partir du
paramètre ``courant`` de ``input.spec``.

**Anisotropie.** Il n'y a pas de type de matériau anisotrope dédié : un milieu
anisotrope (solide ou fluide) se déclare avec son type de base ``S``/``F``, et
l'anisotropie est portée par ``material.spec`` (``deftype`` de type ``*_Aniso`` /
``Cstar*``), exactement comme les autres propriétés variables. L'ancien type ``A``
(fluide anisotrope) est déprécié — il reste accepté comme alias de ``F`` mais n'est
plus émis. Le mailleur lit ``material.spec`` pour graver le bon domaine (fluide
anisotrope) dans le maillage.


.. _pml.input:

Format de pml.input (PML par extrusion)
=======================================

Le fichier ``pml.input`` est **facultatif**. Il permet d'ajouter des couches de
PML à un maillage **importé** (UNV, Abaqus, HDF5) qui n'en contient pas, sans
avoir à régénérer le maillage à la source. C'est le pendant, pour les maillages
externes, des PML « on the fly » définies dans ``mat.dat`` (qui ne concernent
que la grille cartésienne automatique).

En l'absence de fichier ``pml.input``, **aucune PML n'est ajoutée** au maillage
importé (comportement inchangé).

La fonctionnalité est disponible pour les deux mailleurs :

- :program:`mesher` (3D) : côtés ``x-`` ``x+`` ``y-`` ``y+`` ``z-`` ``z+`` ;
- :program:`mesher2D` (2D) : côtés ``x-`` ``x+`` ``z-`` ``z+`` seulement
  (``z`` est l'axe vertical ; ``y-``/``y+`` sont rejetés).

**Sélection de la source des PML (maillages importés, choix 2/3/4).** Le mailleur
lit d'abord ``mater.in`` :

- **si ``mater.in`` déclare déjà des matériaux PML** (type ``P``/``L``), le maillage
  importé contient déjà les éléments PML : le mailleur en **déduit les descripteurs**
  (positions/largeurs et matériau associé) à partir de la géométrie et **ignore
  ``pml.input``** ;
- **sinon**, si ``pml.input`` est présent, le mailleur :

  1. détecte les faces (3D) / arêtes (2D) de bord situées sur le côté demandé ;
  2. les extrude vers l'extérieur du nombre de couches d'éléments indiqué ;
  3. crée les matériaux PML dérivés en copiant les propriétés isotropes du matériau
     adjacent.

Dans les deux cas, **seul ``material.input`` est (ré)écrit**. Les PML y sont des
matériaux isotropes standard (``P``/``L``) ; ``material.spec``, s'il existe, ne
définit que les matériaux intérieurs (aucun bloc PML n'y est ajouté). Une PML est
toujours **isotrope** : elle reprend le Vp/Vs/Rho du matériau de bord, même si
celui-ci est aléatoire ou anisotrope.

Pour un maillage « on the fly » (choix 1), ``pml.input`` est ignoré (les PML
proviennent alors de ``mat.dat``).

Format du fichier
-----------------

Une ligne par côté à traiter, les lignes de commentaire commençant par ``#`` ::

  # Format standardisé (3D ou 2D) :
  # <côté> [épaisseur totale] [nb d'éléments] [ratio/exposant] [loi]
  #
  # Note : En 2D, l'ancien format "<côté> [nb d'éléments] [épaisseur totale]" reste aussi supporté par compatibilité.
  #
  x- 250. 3 1.2 geom       # 3 couches sur 250m d'épaisseur totale, progressant géométriquement avec un ratio de 1.2
  x+ 250. 3 2.0 power      # progressant avec une loi de puissance (power-law) d'exposant 2.0 (finesse proche du domaine)
  y- 250. 3 3.0 linear     # progressant linéairement avec un ratio de taille d'élément final/initial de 3.0
  y+ 250. 3 1.0 geom       # ratio=1.0 ou loi omise => maillage homogène (constant)
  z- 3                    # sans épaisseur (épaisseur auto d'une couche de bord par élément)

- **côté** : ``x-`` ``x+`` ``y-`` ``y+`` ``z-`` ``z+`` (un côté absent = pas de PML).
- **épaisseur totale** : épaisseur totale de la PML sur ce côté.
- **nb d'éléments** : nombre de couches d'éléments PML.
- **ratio** (optionnel, défaut 1.0) : ratio de progression géométrique (geom), d'exposant de puissance (power), ou ratio de taille d'élément final/initial (linear).
- **loi** (optionnelle, défaut "geom") : type de loi de répartition : ``geom`` (géométrique), ``power`` (power-law) ou ``linear`` (linéaire).

En 2D, on peut de plus préciser les paramètres d'atténuation communs à toutes
les PML créées ::

  pmlparams <npow> <Rc> <omegac> <kc>

(Les PML utilisent le type choisi dans ``input.spec`` via la section
``pml_infos { pml_type = PML|CPML|... }``.)

Coins et arêtes
---------------

Les côtés sont traités dans un ordre fixe (x, puis y, puis z). Les coins et les
arêtes sont donc générés **automatiquement** : lorsque la passe en ``y``
rencontre la face ``y-`` d'une colonne PML déjà créée en ``x``, elle produit un
matériau PML de coin combinant les directions (par ex. ``W+S``). Il n'y a rien à
déclarer pour les coins.

Limitations (v2)
----------------

- Les éléments à 8 nœuds (Hexa8 en 3D), 27 nœuds (Hexa27 en 3D), 4 nœuds (Quad4 en 2D) et 8 nœuds (Quad8 en 2D) sont pleinement supportés.
- Le côté choisi doit être un plan aligné sur un axe (il coïncide avec le plan
  de la boîte englobante de ce côté) ; une surface non plane (topographie) ne
  peut pas être extrudée.

Exemples
--------

- 3D : ``SEM/SEM3D/TESTS/NON-REGR/TEST_0011_cube_pml_extrude`` ;
- 2D : ``SEM/SEM2D/TESTS/2Dsquare_pml_extrude`` ;
- mailleur seul (2D) : ``SEM/MESH2D/TESTS/hdf5_pml``.


Format de material.input (version obsolète)
===========================================

**Dans une future version le contenu de ce fichier sera intégré au fichier input.spec**

Le fichier ``material.input`` est créé automatiquement pour le cas avec un maillage automatique.

Pour le cas où le maillage n'est pas automatique, le fichier a par exemple l'allure suivante ::

  22
  P 0380 150 1900 05 05 07 0.000005 0 0
  P 1100 180 1900 05 05 07 0.000005 0 0
  P 1100 180 1900 07 05 07 0.000005 0 0
  P 1100 180 1900 07 07 07 0.000005 0 0
  P 1100 180 1900 05 07 07 0.000005 0 0
  P 1100 180 1900 07 07 07 0.000005 0 0
  P 1100 180 1900 07 05 07 0.000005 0 0
  S 0380 150 1900 05 05 05 0.000005 0 0
  S 1100 180 1900 05 05 05 0.000005 0 0
  P 1100 180 1900 07 05 05 0.000005 0 0
  P 1100 180 1900 07 07 05 0.000005 0 0
  P 1100 180 1900 05 07 05 0.000005 0 0
  P 1100 180 1900 07 07 05 0.000005 0 0
  P 1100 180 1900 07 05 05 0.000005 0 0
  P 0380 150 1900 05 05 07 0.000005 0 0
  P 1100 180 1900 05 05 07 0.000005 0 0
  P 1100 180 1900 07 05 07 0.000005 0 0
  P 1100 180 1900 07 07 07 0.000005 0 0
  P 1100 180 1900 05 07 07 0.000005 0 0
  P 1100 180 1900 07 07 07 0.000005 0 0
  P 1100 180 1900 07 05 07 0.000005 0 0
  S 2000 900 1900 05 05 05 0.000005 0 0
  # PML properties
  # Filtering? npow,Apow1 X+X-Y+Y-Z+Z-
  F 2 10. F F F F T T 0.
  F 2 10. F F F F T T 0.
  F 2 10. T T F F T T 0.
  F 2 10. T T T T T T 0.
  F 2 10. F F T T T T 0.
  F 2 10. T F T T T T 0.
  F 2 10. T F F F T T 0.
  F 2 10. T T F F F F 0.
  F 2 10. T T T T F F 0.
  F 2 10. F F T T F F 0.
  F 2 10. T F T T F F 0.
  F 2 10. T F F F F F 0.
  F 2 10. F F F F T F 0.
  F 2 10. F F F F T F 0.
  F 2 10. T T F F T F 0.
  F 2 10. T T T T T F 0.
  F 2 10. F F T T T F 0.
  F 2 10. T F T T T F 0.
  F 2 10. T F F F T F 0. 

-  Le format du fichier est le suivant :
  
  - la première ligne contient le nombre de milieux décrits
  
  - Une ligne par milieu, contenant :
  
    - le type de milieu (Solide, Fluide, PML solide (P)m PML fluide (L) )
  
    - Les vitesses d'ondes P, et S
  
    - La densité
  
    - L'ordre des éléments en X, Y, Z (Y est ignoré en 2D)
  
    - Un pas de temps (ignoré dans la version actuelle)
  
    - Les attenuations d'ondes P et S
  
  - 2 lignes de commentaires
  
  - Pour chaque milieu de type PML (donc P ou L), une ligne indiquant les directions d'atténuation,
    et le type d'attenuation :
  
    - Un caractère pour le type de PML (filtrante (T), ou standard (F))
  
    - paramètres n et A pour les PML filtrantes
  
    - 3 couples de deux drapeaux T ou F (pour True False) indiquant si la PML atténue dans
      les directions X, Y et Z respectivement (premier flag du couple) et dans le sens positif (T)
      ou négatif de l'axe.
  
    - La fréquence de coupure en cas de PML filtrante



Format de material.input (nouvelle version)
===========================================

**Format commun à SEM2D et SEM3D** (depuis l'unification du format 2D avec le
3D). En 2D, l'axe vertical est ``z`` : les champs ``posY``/``widthY`` du bloc PML
valent toujours 0, et les directions d'atténuation (gauche/droite, haut/bas) sont
déduites du signe des largeurs ``widthX``/``widthZ``.

Le format de ``material.input`` a été modifié.
Voici un exemple de sa nouvelle présentation, suivi d'explications sur les paramètres qui interviennent : ::

    21
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    S 630.000000 250.000000 1800.000000 0.000000 0.000000
    S 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    P 630.000000 250.000000 1800.000000 0.000000 0.000000
    # PML properties
    # npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat
    2 10. -2.450000 0.000000 0.000000 0.000000 6.000000 -5.000000 7
    2 10. -5.000000 0.000000 -5.000000 0.000000 6.000000 -5.000000 8
    2 10. -5.000000 -5.000000 -5.000000 0.000000 6.000000 -5.000000 8
    2 10. -5.000000 -5.000000 -5.000000 -5.000000 6.000000 -5.000000 8
    2 10. -5.000000 0.000000 -5.000000 -5.000000 6.000000 -5.000000 8
    2 10. 5.000000 5.000000 -5.000000 -5.000000 6.000000 -5.000000 8
    2 10. 5.000000 5.000000 -5.000000 0.000000 6.000000 -5.000000 8
    2 10. -5.000000 -5.000000 -5.000000 0.000000 6.000000 0.000000 8
    2 10. -5.000000 -5.000000 -5.000000 -5.000000 6.000000 0.000000 8
    2 10. -5.000000 0.000000 -5.000000 -5.000000 6.000000 0.000000 8
    2 10. 5.000000 5.000000 -5.000000 -5.000000 6.000000 0.000000 8
    2 10. 5.000000 5.000000 -5.000000 0.000000 6.000000 0.000000 8
    2 10. -2.450000 0.000000 0.000000 0.000000 26.000000 5.000000 7
    2 10. -5.000000 0.000000 0.000000 0.000000 26.000000 5.000000 8
    2 10. -5.000000 -5.000000 0.000000 0.000000 26.000000 5.000000 8
    2 10. -5.000000 -5.000000 -5.000000 -5.000000 26.000000 5.000000 8
    2 10. -5.000000 0.000000 -5.000000 -5.000000 26.000000 5.000000 8
    2 10. 5.000000 5.000000 -5.000000 -5.000000 26.000000 5.000000 8
    2 10. 5.000000 5.000000 -5.000000 0.000000 26.000000 5.000000 8

Plus précisément, ce fichier doit contenir :
  
- Le nombre de milieux décrits sur la première ligne ;
  
- Une ligne par milieu, contenant :
  
  - Le type de milieu :
        
        - Solide : noté S
        
        - Fluide : noté F
        
        - PML solide : noté P
        
        - PML fluide : noté L
  
  - Les vitesses d'ondes P et S (en m/s)
  
  - La densité (en kg/m^3)
  
  - Les coefficients attenuations d'ondes P et S.
  
- 2 lignes de commentaires
  
- Pour chaque milieu de type PML (donc P ou L), une ligne permettant de paramétrer les directions d'atténuation :
  
  - Les paramètres n et A pour les PML filtrantes ;
  
  - Les coordonnées (posX, posY, posZ) d'un point à la frontière entre le PML et la matériau solide/liquide adjacent ;
  
  - La taille de l'extrusion suivant la direction de l'espace (widthX, widthY, widthZ). Elle peut être positive ou négative suivant le sens donné aux axes (X,Y,Z). Une inversion du signe occasionnera une amplification au lieu d'une atténuation !
  
  - Le numéro attributé au matériau solide/liquide adjacent au PML. Attention ! la numérotation commence à 0. Ainsi les deux matériaux solides décrits dans l'exemple ci-dessus sont bien numérotés 7 et 8 (et non 8 et 9).

**Remarque :** pour l'extrusion, trois types de figure sont à distinguer, suivant la nature de la surface de contact PML-solide/liquide :

- Surface de contact : l'extrusion se fait uniquement dans la direction normale à la surface de contact. On mettra 0 dans les deux autres directions.
   
- Arête en contact : l'extrusion se fait dans les deux directions normales à l'arête uniquement.
   
- Coin en contact : si seul un point fait le contact entre le PML et le milieu adjacent, l'extrusion se fait dans les trois dimensions de l'espace.
   
  
   
Format de material.spec
=======================

.. _material.spec:

Le fichier ``material.spec`` est facultatif [#]_ et peut contenir des définitions
concernant les matériaux utilisés lors du calcul.

.. [#] Dans une future version, ``material.input`` va disparaître et ``material.spec``
   deviendra obligatoire.


Syntaxe du fichier
------------------

La syntaxe générale est la même que celle du fichier ``input.spec``. Une seule
section est valide pour l'instant. Voici un exemple ::

  material 0 {
     domain = solid;
     deftype = Kappa_Mu_Rho;
     spacedef = file;
     filename0 = "mat/h5/Mat_0_Kappa.h5";
     filename1 = "mat/h5/Mat_0_Mu.h5";
     filename2 = "mat/h5/Mat_0_Density.h5";
  };

  material 1 { copy = 0; };

Chaque section ``material`` est suivie du numéro du matériau concerné.

Le contenu de la section est composé de différentes variables :

=============  ===================================================================
Nom            Description
=============  ===================================================================
domain         Domaine de calcul associé (solid|fluid|solidpml|fluidpml) [#]_
deftype        Indique quelles sont les variables utilisées pour la définition
spacedef       Indique si les propriétés sont variables ou constante spatialement
filename       Nom d'un fichier contenant toutes les variables
filename0      Nom d'un fichier contenant la première variable
filename1      Nom d'un fichier contenant la deuxième variable
filename2      Nom d'un fichier contenant la troisième variable
copy           Numéro d'un matériau dont on copie la définition [#]_
Vp             Si ``spacedef=constant`` : vitesse d'onde P
Vs             Si ``spacedef=constant`` : vitesse d'onde S
Rho            Si ``spacedef=constant`` : densité

=============  ===================================================================

Description des mot-clefs de ``deftype`` :

- ``Kappa_Mu_Rho`` : Variable 0 : kappa, Variable 1 : Mu, Variable 2 : Rho
- ``Lambda_Mu_Rho`` : ...
- ``Vp_Vs_Rho``
- ``E_nu_Rho`` : Module d'Young, coefficient de Poisson, Densité
- ``Hooke_Rho`` : Cijkl, Rho.
- ``CStar``: input pour le code homofft (solide anisotrope, 21 composantes + Rho).
- ``Fluid_Aniso``: tenseur de module volumique anisotrope :math:`K_{ij}` + densité (7 propriétés),
  lu depuis un fichier HDF5 avec des groupes ``K11``, ``K22``, ``K33``, ``K12``, ``K13``,
  ``K23``, ``Rho``. Voir :ref:`anisotropic_fluid`.
- ``Cstar_Fluid``: format binaire Cstar de homofft restreint au cas acoustique (Nd=3),
  7 composantes par point (K11, K12, K13, K22, K23, K33, Rho).
  Voir :ref:`anisotropic_fluid`.

Description des mot-clefs de ``spacedef`` :

- ``constant`` : les valeurs sont précisées dans la suite de la section ``material``
- ``file`` : les valeurs sont données dans un ou plusieurs fichiers.

.. [#] Evolution prévue : on conservera solid|fluid, la qualification de PML
   sera gérée dans le maillage directement

.. [#] La principale utilitée actuellement est de donner les même propriétés aux
   matériau PML. Une fois la modification précédente implémentée, cela pourra servir
   à distinguer deux zones de matériaux identiques (pour les sorties par exemple)


Format des fichiers materiaux
-----------------------------

Les fichiers spécifiés par ``filename`` donnent les propriétés matériaux sur une
grille régulière. Le format du fichier est HDF5 avec une structure imposée.

Le code teste deux possibilités pour la structure. Cela permet de stocker les
variables dans un seul fichier ou dans un fichier par variable.

Dans le premier cas, le fichier doit contenir un groupe par variable, le nom
du groupe devant être le nom de la variable (soit pour l'instant : ``Vp``, ``Vs``, ``Rho``,
``E``, ``Nu``, ``Lambda``, ``Mu``, ``Kappa``).

Dans le second cas (un fichier par variable), deux options sont possibles : Soit la
variable est stockée dans un groupe portant son nom (comme dans le premier cas), soit
elle est stockée à la racine du fichier.

La description d'une variable (à la racine ou dans un groupe) est la même : le fichier ou groupe
doit contenir :

- ``xMinGlob`` : un attribut de taille 3 contenant les coordonnées minimum de la grille

- ``xMaxGlob`` : un attribut de taille 3 contenant les coordonnées maximum de la grille

- ``samples`` : un dataset comportant 3 dimensions, la taille de chaque dimension étant de 2 minimum.

Le dataset est stocké dans l'ordre Fortran (nx,ny,nz) donc apparaîtra comme ayant une taille de (nz,ny,nx) depuis un code C (ou avec ``h5dump``).

La valeur au point (0,0,0) du tableau ``samples`` correspond au coordonnées spatiales (xMinGlob(0), xMinGlob(1), xMinGlob(2)).

La valeur au point (nx,ny,nz) (ordre fortran) correspond au coordonnées spatiales (xMaxGlob(0), xMaxGlob(1), xMaxGlob(2)).

Les trois variables ne sont pas nécessairement définies sur la même grille.

