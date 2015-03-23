.. -*- coding: utf-8; mode: rst -*-

============
Introduction
============


SEM est un code de propagation d'onde sismique, en 2D et 3D, fondé sur la méthode
des éléments spectraux ([PRI94]_, [FAC97]_, [KOM98]_, [SER98]_, [KOM99]_).

Cette version prend en compte la propagation dans des milieux
hétérogènes à géométrie complexe (topographie de surface, interfaces),
de type :

- solide : élastique ou viscoélastique,  isotrope ou anisotrope, aléatoire ou non,

- ou fluide.

Les conditions d'absorption en bord de domaine utilisent la méthode
*Perfectly Matched Layer (PML)* ([BER94]_, [FES05]_).

Un outil de prétraitement permet de convertir et partitionner (avec
:program:`Metis`) des maillages au format :program:`Abaqus` ou
:program:`Ideas` (UNV), ou de créer des maillages cartésiens
structurés.

Le code est écrit majoritairement en Fortran 90. Il utilise les librairies `BLAS <http://www.netlib.org/>`_, `MPI <http://www.openmpi.org>`_, et `HDF5 <http://www.hdfgroup.org>`_.

.. [BER94] Bérenger, J.-P. (1994). A perfectly matched layer for the absorption of electromagnetics waves, *J. Comp. Phys. 114*, 185-200.

.. [FAC97] Faccioli, E., R. Maggio, R. Paolucci, and A. Quarteroni (1997). 2D and 3D elastic wave propagation by a pseudo-spectral domain decomposition method. *J. Seismol. 1*, 237-251.

.. [FES05] Festa, G. and J.-P. Vilotte (2005). The Newmark scheme as velocity-stress time-staggering: an efficient PML implementation for spectral element simulations of elastodynamics. *Geophys. J. Int. 161(3)*, 789-812.

.. [KOM98] Komatitsch, D. and J.-P. Vilotte (1998). The spectral element method: an efficient tool to simulate the seismic response of 2D and 3D geological structures. *Bull. Seismol. Soc. AM. 88(2)*, 368-392.

.. [KOM99] Komatitsch, D. and J. Tromp (1999). Introduction to the spectral-element method for three-dimensional seismic wave propagation. *Geophys. J. Int. 139(3)*, 806-822.

.. [PRI94] Priolo, E., J. M. Carcione, and G. Seriani (1994). Numerical simulation of interface waves by high-order spectral modeling techniques. *J. Acoust. Soc. AM. 95(2)*, 681-693.

.. [SER98] Seriani, G. (1998). 3D large-scale wave propagation modeling by a spectral element method on a Cray T3E multiprocessor. *Computer Methods in Applied Mechanics and Engineering 164(1)*, 235-247.




Notations
---------

Exemple des notations utilisées dans ce document :

- Un programme ou logiciel : :program:`CMake`

- Une commande Unix (``cmake``), un nom de fichier (``/usr/bin/python``) ou un texte à saisir ("saisissez le mot-clef ``false``").

- Un terme en anglais : *remote control*

- Une variable d'environnement : :envvar:`HDF5_ROOT`

Parfois le programme et sa commande Unix ont le même nom, on essaira de faire la distinction, par exemple :

"Pour configurer un programme utilisant :program:`CMake` il faut taper la commande ``cmake``."


Quoi de neuf
------------

.. toctree::

   NEWS.rst

Evolutions futures
------------------

Certaines fonctionnalités sont prévues (voire déjà disponibles dans le code) mais
n'ont pas encore été finalisées, intégrées ou correctement testées :

- Description de gradient de propriétés dans les matériaux. Le code de la version CEA
  a été intégré, mais la description des matériaux dans le fichier de configuration
  n'a pas encore été effectuée.

  La nouvelle description des gradients et le nouveau format du fichier matériaux
  seront développés dans une future version.

- Description des conditions de Neumann. Le code existe, il n'a pas été testé. Il sera intégré
  dans le fichier de configuration au nouveau format dans une prochaine version.

- Anisotropie : le code pour gérer des matériaux anisotropes existe,
  mais il n'y a rien dans la syntaxe actuelle du fichier de
  description des matériaux qui permette de définir un milieu
  anisotrope. Là encore, cela sera intégré dans la prochaine version
  lors de la refonte du fichier de description des matériaux.



Notes importantes
-----------------

Le code source est versionné avec :program:`Git` et livré dans une archive contenant :

- SEM version 3D

- SEM version 2D

- MESH : un outil de préparation de maillages 3D pour :program:`SEM3D` (l'équivalent
  2D sera intégré dans une prochaine version).

- La librairie :program:`HDF5` est devenue une dépendance obligatoire (
  `www.hdfgroup.org <http://www.hdfgroup.org>`_ ).

  Cette librairie permet le stockage efficace de gros volume de
  données. Son utilisation permet le post-traitement immédiat des
  snapshot avec :program:`Paraview` ou :program:`Ensight`. Les données produites sont
  également lisibles facilement avec :program:`Matlab` et :program:`Python`.

- Le schéma en temps a été simplifié (Les paramètres beta/gamma de
  l'algorithme de Newmark ne sont plus modifiables).

  Ils pourront être réintroduits une fois réglé le problème de
  synchronisation avec les forces de couplage externes.

- Bien que les deux méthodes continuent de coexister, le calcul des
  forces utilisant le tableau ``Acoeff`` a été désactivé dans cette
  version. Le code est plus lisible mais moins rapide.

  On étudiera comment obtenir le meilleur des deux méthodes dans une
  prochaine version.


