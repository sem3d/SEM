.. -*- mode:rst; coding: utf-8 -*-

Maillage d'une topographie
==========================


Génération du maillage
----------------------

Pour générer ce cas on va utiliser un jeu d'outils externes à SEM : *meshtools*.

Les étapes de construction sont les suivantes :

- Sélection d'un ou plusieurs fichiers de topographie (format SRTM par exemple) (*utilisateur*)

- Conversion/concaténation de la topographie en un format compact intermédiaire (*mt_import*)

- Création d'une grille cartésienne dans la projection souhaitée (*mt_grid*)

- (optionnel) Création de grilles supplémentaires pour mailler des couches en profondeur épousant la topographie
  de surface (*utilisateur*)

- Génération du maillage et du fichier matériau associé (*mt_topo*)

- Partition du maillage (*mesher*)


Nous allons traiter un exemple de génération de maillage à partir d'un fichier srtm ::

  # On décompresse le fichier srtm
  $ unzip srtm_56_01.zip
  # On convertit le fichier au format hdf5 (lat/lon)
  $ mt_import -s topo_srtm.h5 srtm_56_01.tif
  # On projete une grille de 30x30 mailles de 1000x1000 m de cote d'origine 58N 96E dans la projection aeqd
  $ mt_grid --vx=1000,0 --vy=0,1000 -g 30,30 -p "+proj=aeqd +lat_0=58.0 +lon_0=96.0" -n surf topo_srtm.h5 grid.h5
  $ mt_grid --vx=500,0 --vy=0,500 -g 300,300 -p "+proj=aeqd +lat_0=58.0 +lon_0=96.0" -n surf topo.h5 grid.h5
  # Le fichier contenant la grille est utilise pour creer un maillage
  $ mt_topo --npml=1 --profile=mesh.profile --mat=input_material.dat grid.h5 mesh_sem.h5
  # on renomme le fichier materiau (pour l'outil mesher)
  $ cp mesh_sem.h5.mat material.input
  $ mesher
  256
  0
  4
  1
  mesh_sem.h5
  $ mkdir sem
  $ mv mesh4spec.0* sem/
  # Lancement du cas sem
  $ mpirun -n 256 sem3d.exe

Modification de l'association des matériaux
-------------------------------------------

L'outil ``mt_topo`` via le fichier de profil vertical (option ``--profile``) applique une description
de milieu homogène par couche de mailles (pas de variation en X et Y).

On peut cependant aller plus loin et modifier le maillage généré avec quelques lignes de script python ::

  $ python
  # import des fonctions numpy
  >>> from numpy import *
  # Import du module de lecture de fichier HDF5
  >>> import h5py
  # Ouverture du fichier
  >>> fmesh = h5py.File("mesh_sem.h5","r+")
  # On lit les coordonnees des noeuds (taill Np x 3)
  >>> nodes = fmesh["/Nodes"][...]
  # On charge les proprietes materiau (taille Nel)
  >>> mat = fmesh["/Mat"][...]
  # On charge la description des elements Nel x 8
  >>> elem = fmesh["/Elements"][...]
  # On calcule le centre de chaque element nodes[elem,:] est un tableau
  # de taille Nel x 8 x 3, on fait la moyenne des coordonnees sur l'axe du milieu
  >>> ctr = nodes[elem,:].sum(axis=1)/8.
  # on applique un nouveau materiau sur la zone d'interet :
  >>> z1 = logical_and( ctr[:,0] > 5000, ctr[:,0] < 10000. )
  >>> z2 = logical_and( ctr[:,1] > 2000, ctr[:,1] < 4000. )
  >>> z3 = ctr[:,2] > -5000
  # Un tableau de booléen de taille Nel tq les valeurs true correspondent aux
  # élements de centre 5000<X<10000 , 2000<Y<4000, Z>-5000
  >>> zone = logical_and(z1, logical_and(z2, z3))
  # On change le materiau associé à cette zone
  >>> mat[zone] = 2
  # On récrit le nouveau champ matériau
  >>> fmesh["/Mat"] = mat
  # Fin
  >>> fmesh.close()
