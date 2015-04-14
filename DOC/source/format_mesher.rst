.. -*- coding: utf-8 -*-

============================
Format des fichiers matériau
============================

.. _material.input: 

Format de mat.dat
=================

Ce fichier est seulement nécessaire dans le cas où le maillage est realisé par le mailleur automatique.
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

Format de mater.in
==================

Le fichier ``mater.in`` décrit combien de matériaux sont utilisés dans le modèle :: 

  1
  S  6300.00  2500.00   2800. 5   5    5  0.000005 630. 250.

 `1` est le nombre de matériaux dans le modèle.

La deuxième ligne décrit le type de matériau (``S`` matériau solide et
``F`` matériau fluide).Pour chaque matériau, on déclare
successivement, la vitesse de propagation de l'onde de pression,
vitesse de l'onde de cisaillement, la densité du matériau, le nombre
de gll dans les trois directions( ``x``, ``y`` et ``z``), le pas de
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
  
  - la première ligne contient le nomnbre de milieux décrits
  
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
  
    - 3 couples de deux drapeaux T ou F (pour True False) indiquant si la PML attenue dans
      les directions X, Y et Z respectivement (premier flag du couple) et dans le sens positif (T)
      ou négatif de l'axe.
  
    - La fréquence de coupure en cas de PML filtrante
