.. -*- coding: utf-8 -*-

============================
Format des fichiers matériau
============================


Format de mat.dat
=================

Cette fichier est juste pour le cas oú le maillage est realise pour le ``mesher`` 
automatique.
Le fichier ``mat.dat`` doit contenir (les commentaires, après le *#*
sont facultatifs) ::

  -100.    # xmin
  500.     # xmax
  50.      # xstep
  -100.    # ymin
  500.     # ymax
  50.      # ystep
  500.     # zmax
  1        # nb. of layers
  600. 12  # upper layer: thickness and nb of steps
  1  # PMLs? 0: no, 1: yes
  1 1  # PMLs on top? at the bottom? (0: no, 1: yes)
  5   # nb of GLL nodes in the PML
  1   # 8 or 27 control points for elements (1 or 2)

- Choix de 8 noeuds par maille : 1 (Les mailles quadratiques à 27
  noeuds sont en développement)

Format de mater.in
==================

Le fichier ``mater.in`` décrit combien de matériels ont dans le model :: 

  1
  S  6300.00  2500.00   2800. 5   5    5  0.000005 630. 250.

 `1` est le nombre de matériaux dans le model.

La deuxième ligne décrit le type de matériau (``S`` matériau solide et
``F`` matériau fluide).Pour chaque matériau, on déclare
successivement, la vitesse de propagation de l'onde de pression,
vitesse de l'onde de cisaillement, la densité du matériau, le numbre
de ggl dans les trois directions( ``x``, ``y`` et ``z``), le pas de
temps (ignoré dans la version actuelle) et les atténuations des ondes
P et S.


Format de material.input
========================

**Dans une future version le contenu de ce fichier sera intégré au fichier input.spec**

Le fichier ``material.input`` est créé automatiquement pour le cas avec un maillage automatique.

Pour le cas où le maillage n'est pas automatique le ficher doit contenir ::

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

Les lettres ``T`` et ``F`` (True et False) sont utilisées pour définir les directions de l'atténuation de la PML. Trois couples de deux 
drapeaux T ou F indiquant si la PML attenue dans les directions X, Y et Z respectivement 
(premier flag du couple) et dans le sens positif (T) ou négatif de l'axe.
