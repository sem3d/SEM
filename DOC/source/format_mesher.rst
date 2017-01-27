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
  5        # nb of GLL nodes in the PML
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
  5        # nb of GLL nodes in the PML
  1        # 8 or 27 control points for elements (1 or 2)
  

Format de mater.in
==================

Le fichier ``mater.in`` décrit combien de matériaux sont utilisés dans le modèle :: 

  1
  S  6300.00  2500.00   2800. 5   5    5  0.000005 630. 250.

 `1` est le nombre de matériaux dans le modèle.

La deuxième ligne décrit le type de matériau (``S`` matériau solide et
``F`` matériau fluide). Pour chaque matériau, on déclare
successivement, la vitesse de propagation de l'onde de pression,
vitesse de l'onde de cisaillement, la densité du matériau, le nombre
de GLLs dans les trois directions( ``x``, ``y`` et ``z``), le pas de
temps (ignoré dans la version actuelle) et les paramètres :math:`Q_\kappa` et :math:`Q_\mu`
pour l'atténuation des ondes P et S.


Format de material.input
========================

**Dans une future version le contenu de ce fichier sera intégré au fichier input.spec**

Le fichier ``material.input`` est créé automatiquement pour le cas avec un maillage automatique.

Pour le cas où le maillage n'est pas automatique, le ficher doit contenir ::

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



Format de material.input V2
===========================


::

   Nmat
   Type  Pspeed  Sspeed Density  NGLL  Qp  Qmu
   # Two lines skipped after Nmat material lines
   # The second line
   npow  Apow  posX  widthX  posY widthY posZ widthZ  mat


Description des paramètres : comme le format V1 sauf :

- posX, widthX : position du début de la couche PML sur l'axe X. Si widthX == 0, alors
  pas de PML en X. Le signe de widthX donne sa direction : la pml
  se trouve entre (posX) et (posX+widthX)

   
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
- ``Hooke_Rho`` : (Non-implémentée) Cijkl, Rho.
  
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

