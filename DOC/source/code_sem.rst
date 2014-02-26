.. -*- coding: utf-8 -*-

.. _code_sem:

=======================
Description du code SEM
=======================


SEM est un code éléments finis de simulation de la propagation d'ondes
dans un milieu solide ou fluide (ou mixte).

Dans les deux milieux, l'équation d'onde est résolue par la méthode
des éléments spectraux : les champs (déplacement, vitesse, ...) sont
évalués sur des fonctions de bases de chaque élément du maillage. Dans
le cas de SEM, ces fonctions sont des polynômes de Lagrange, dont les
noeuds sont placés aux points de Gauss-Lobatto-Legendre.

C'est le placement particulier de ces points qui confère à la méthode
une convergence spatiale rapide en ordre de l'élément.

