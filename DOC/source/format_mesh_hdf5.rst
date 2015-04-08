.. -*- coding: utf-8 -*-

.. _Format HDF5:

===================================
Format du maillage HDF5 pour mesher
===================================

Le partitionneur de SEM accepte un format spécifique pour décrire un maillage en entrée.

La structure sous-jacente est un fichier au format HDF5. Ce fichier doit contenir trois
*dataset* spécifiques :

- ``Elements`` : un tableau de NE x 8 entiers de 0 à (NN-1) faisant
  référence aux noeuds. Les noeuds doivent respecter un ordre spécifique : si on visualise un cube,
  les noeuds 1,2,3,4 sont les noeuds des sommets du carré à la base du cube, dans l'ordre trigonométrique.
  Les noeuds 5,6,7,8 sont les noeuds des sommets du carré au sommet du cube, dans l'ordre trigonométrique, le noeud 5
  se situant directement au dessus du noeud 1.

  Ainsi le noeud 5 est au-dessus du noeud 1, le 6 au dessus du 2, etc...

- ``Nodes`` ; un tableau de ( NN x 3 ) rééls : les coordonnées des
  noeuds

- ``Mat`` : un tableau de NE x 2 entiers. Le premier entier détermine le numéro du matériau à
  associer à chaque maille. Le second est pour l'instant inutilisé, il est réservé pour déterminer
  dans le futur si une maille est une maille PML est dans quelle direction elle absorbe.


Voici ci dessous un exemple de script Python qui permet de créer un cube ::

   # import des fonctions numpy
   import numpy as np
   # Import du module de lecture/écriture de fichier HDF5
   import h5py
   # Ouverture du fichier
   fmesh = h5py.File("mesh_sem.h5","w")
   # On cree un tableau contenant les sommets du cube
   nodes = np.array([
        [0,0,0],
	[1,0,0],
	[1,1,0],
	[0,1,0],
        [0,0,1],
	[1,0,1],
	[1,1,1],
	[0,1,1]], float)
   fmesh["/Nodes"] = nodes
   # On cree le tableau contenant les elements (ici 1 seul)
   elems = np.zeros( (1,8), int)
   elems[0,:] = (0,1,2,3,4,5,6,7)
   fmesh["/Elements"] = elems
   # On definit les proprietes materiau (taille Nel)
   fmesh["/Mat"] = np.array( [[0, 0]], int)
   # Fermeture du fichier
   fmesh.close()
