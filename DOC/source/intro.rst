.. -*- coding: utf-8; mode: rst -*-

============
Introduction
============


SEM est un code de propagation d'onde sismique, en 2D et 3D, fondé sur la méthode
des éléments spectraux ([PRI94]_, [FAC97]_, [KOM98]_, [SER98]_, [KOM99]_).

Cette version prend en compte la propagation dans des milieux hétérogènes, à géométrie complexe (topographie
de surface, interfaces).

Un outil de prétraitement, permet de convertir et partitionner (avec :program:`Metis`) des maillages au format :program:`Abaqus`
ou :program:`Ideas` (UNV).

Des sources ponctuelles ou planaires peuvent être introduites. Les conditions d'absorption en bord de domaine
utilisent la méthode *Perfectly Matched Layer (PML)* décrite dans [FES05]_.

Le code est écrit majoritairement en Fortran 90. Il utilise les librairies `BLAS <http://www.netlib.org/>`_, et `MPI <http://www.openmpi.org>`_.

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

- Une commande Unix : ``cmake``

- Un terme en anglais : *remote control*

- Une variable d'environnement : :envvar:`HDF5_ROOT`