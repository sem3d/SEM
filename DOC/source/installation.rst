.. -*- coding: utf-8 -*-

.. _installation:

==========================
Installation du code SEM3D
==========================


Prérequis
=========

Le code :program:`SEM3D` nécessite deux (ou trois) outils externes pour sa compilation :

- :program:`CMake` : `www.cmake.org <http://www.cmake.org>`_ est un générateur de Makefile (comme
  autoconf/autotools).

- :program:`HDF5` : `www.hdfgroup.org <http://www.hdfgroup.org>`_ est une librairie permettant la
  gestion de fichier de données de grande taille.

  Si :program:`HDF5` est installé dans un chemin non-standard, utilisez la
  variable d'environnement :envvar:`HDF5_ROOT` pour indiquer aux scripts
  :program:`CMake` de :program:`SEM` où elle se
  trouve. :envvar:`HDF5_ROOT` doit pointer sur le préfixe
  d'installation, c'est à dire le chemin contenant :file:`bin/`,
  :file:`include/`, etc...

  Les scripts de configurations :program:`CMake` de :program:`SEM` utilise la commande
  ``h5cc -show`` pour détecter le paramètrage de la librairie :program:`HDF5`.

- une librairie :program:`MPI` (:program:`OpenMPI` recommandée, mais :program:`Mpich` ou :program:`Intel MPI` doivent être compatibles).


Compilation
===========

La préparation de la compilation s'effectue avec :program:`CMake`
(commandes: ``ccmake`` ou ``cmake``), ensuite la compilation elle-même
s'effectue avec la commande ``make``.

Organisation des répertoires
----------------------------

Nous recommandons l'organisation suivante (les noms des répertoires sont à titre indicatif) :

- :file:`sem_src` : répertoire contenant les sources,

- :file:`sem_build` : répertoire contenant les binaires (en dehors du répertoire source),

- :file:`sem_debug` : (facultatif) un second répertoire pour une compilation en mode debuggage de SEM.

Préparation de la compilation
-----------------------------

La préparation se fait à l'aide de la commande suivante ::

  $ cd sem_build
  $ ccmake ../sem_src

``ccmake`` est une commande interactive de :program:`CMake` permettant de
paramétrer la compilation. Le paramétrage s'effectue en deux étapes :

- La première étape (configuration) :program:`CMake` va rechercher les
  librairies nécessaires (HDF5 et OpenMPI)

  Lors de cette étape on peut changer le contenu des variables
  affichées (chemin des librairies, nom du compilateur, options de
  compilation).

  **Important** : C'est à cet endroit qu'il faut préciser le mode de
  compilation par la variable : ``CMAKE_BUILD_TYPE``. On peut saisir :
  ``DEBUG``, ``RELEASE``, ou ``RELWITHDEBINFO``. Si on ne saisit
  rien, :program:`SEM` sera compilé avec les options par défaut du compilateur
  (sans optimisation et sans debuggage avec :program:`gcc`, optimisé sans
  debuggage avec :program:`ifort`).

  Lorsqu'on change des variables, il faut reconfigurer (touche ``c``).

- La seconde étape, la génération des fichiers ``Makefile`` ne peut se faire que si
  l'option ``g`` (*generate and exit*) apparait dans
  l'interface. Cette option n'apparait que si la dernière étape de
  configuration n'a pas modifié de variables.

  En effet, il se peut qu'une reconfiguration change d'autres
  variables (lorsqu'on change le compilateur par exemple), il faut
  alors lancer la configuration une seconde fois.

  Lorsque l'étape de configuration ne modifie aucune variable, on peut
  générer les Makefile (touche ``g``).

Compilation
-----------

Une fois la génération terminée, la compilation se fait simplement par la commande ``make``.

Quelques variantes :

- ``make help`` : affiche toutes les cibles possibles.

- ``make -j N`` : compile en parallèle avec N processus (on peut en
  général utiliser N=nombre de processeurs + 1 ou 2).

- ``make -j N -k`` : compile le plus possible sur N processeur. Ne
  s'arrête pas à la première erreur de compilation.

- ``make VERBOSE=1`` : affiche les lignes de commandes exécutées lors de la compilation.


La compilation produit plusieurs exécutables :

- ``build_src/SEM2D/sem2d.exe`` : Code SEM2D.

- ``build_src/SEM3D/sem3d.exe`` : Code :program:`SEM3D`.

- ``build_src/MESH/mesher`` : Outil pour le partitionnement des maillages et la génération au format :program:`HDF5`.


Exécution
---------

Exécution de :program:`SEM3D` monoprocesseur : ::

  $ cd rep_du_cas
  $ ${chemin_build}/SEM3D/sem3d.exe

Exécution de :program:`SEM3D` en MPI : ::

  $ cd rep_du_cas
  $ mpirun -n 4 ${chemin_build}/SEM3D/sem3d.exe

Lancement du générateur de maillage : ::

  $ cd rep_du_cas
  $ ${chemin_build}/MESH/mesher

Ou en mode automatique avec les saisies clavier enregistrées dans le fichier ``mesh.input`` (c'est le cas des cas tests présent avec les sources de SEM) : ::

  $ cd rep_du_cas
  $ ${chemin_build}/MESH/mesher < mesh.input
  

