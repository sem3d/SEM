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
  

Résolutions des problèmes de compilation
----------------------------------------

Dans l'ordre :

1. Lire le message d'erreur

2. Déterminer si il s'agit d'une erreur de compilation ou d'une erreur d'édition de lien

3. **Relire le message d'erreur** et **tout** le message...

4. Regarder ci-dessous si c'est un problème courant


Plusieurs problèmes peuvent survenir lors de la compilation, et/ou l'édition de lien de SEM.

Pour les résoudre il faut avant tout comprendre le processus de compilation :

- Chaque fichier source (``.f``, ``.c``, ``.f90``) est transformé par
  le *compilateur* en un fichier binaire (``.o``).

- En supposant que la version que vous compilez a déjà été compilée
  par ailleurs, les erreurs qui peuvent survenir lors de la
  compilation sont :

  - Un compilateur non testé : Fortran est un langage très mal
    normalisé, ``gfortran`` est souvant plus strict que ``ifort``,
    certaines formulation vont compiler avec l'un et pas avec l'aure.

    Exemple notoire : ifort accepte une structure ``allocatable``
    comme membre d'une autre structure alors que gfortran va exiger un
    ``pointer``

  - L'unité de compilation (fichier .f90 par exemple) utilise un
    module externe non reconnu.

    Exemple classique : gfortran ne peut pas charger le module mpi ou hdf5.

    plusieurs cas :

    - Le module n'est simplement pas trouvé :

      - il faut d'abord trouver le module (``mpi.mod`` ou ``hdf5.mod``)
      
      - puis faire en sorte que le compilateur le trouve : il faut
        ajouter l'option ``-I/chemin/vers/module`` dans une des
        variables ``*_FLAGS`` de cmake. On peut vérifier ce que
        ``cmake`` passe au compilateur avec ``make VERBOSE=1``

    - Le module est produit par un autre compilateur :

      - gfortran ne peut pas utiliser un module compilé avec ifort et vice-versa.

        Il ne peut pas utiliser non plus la librairie produite,
        autrement dit un appel à une fonction externe (dans le style
        Fortran 77) va compiler mais produire des erreurs lors de
        l'exécution.

      - Plus génant, gfortran est incapable d'utiliser un module
        produit par une version majeur différente : on ne peut pas
        compiler SEM avec gfortran 4.8 et lui faire utiliser une
        librairie compilée avec gfortran 4.7

Après la compilation vient l'édition de lien, c'est le moment ou l'on
assemble les fichiers ``.o`` pour en faire un exécutable. Cela
consiste principalement à relier les appels de fonctions externes à
une unité avec leurs définitions dans une autre unité ou dans une
librairie.

Il y a encore plusieurs erreurs classiques :

- La librairie n'est pas trouvée :

  Il faut inclure la librairie dans la compilation. Dans ``cmake`` ce
  sont les variables ``*_LDFLAGS`` ou ``*LIBRARIES`` qui contrôle
  cette partie de la procédure. On peut ajouter le chemin complet
  d'une librairie, ou les options ``-L/chemin -lnom_de_lib``.

  Si le linker indique qu'il ne trouve pas une librairie, c'est que
  celle-ci lui à été désignée : donc soit une option ``-lnom_de_lib``
  existe mais aucun fichier ``libnom_de_lib.so`` n'est présent dans
  les chemins fournis au linker, soit une librairie utilisée indique
  qu'elle dépend d'une autre librairie introuvable.

- La seconde erreur possible est que le linker ne peut pas résoudre un
  symbole. C'est à dire que quelque part dans un fichier ``.o`` ou
  dans une librairie, une fonction est appelée, mais la définition de
  cette fonction n'est dans aucun autre fichier ``.o`` ou librairie.

  - Cela peut venir du code : par exemple lorsqu'on compile SEM en
    monoprocesseur, on utilise une "fausse" librairie MPI. Cette
    version étant moins testée, il se peut qu'un développeur ait
    utilisé une fonction MPI non encore émulée. (c'est arrivé dans le
    passé, mais le nombre de fonction MPI étant limité, on va finir
    par les avoir toutes émulées).

  - Le plus souvent, cela vient d'une librairie qui dépend d'une autre librairie, qui ne spécifie pas
    ses dépendances (car ce n'est pas obligatoire pour une librairie).

    C'est le cas avec les librairies de support du compilateur intel.

    Par exemple: je compile hdf5 avec icc 10. qui utilise des
    fonctions provenant de la librairie de support du compilateur
    10.0. Et plus tard je compile *et* je link SEM avec
    icc/ifort 11. Entre temps la fonction utilisée par intel 10. à
    disparu est n'est plus dans la librairie de support de
    Intel 11. Donc à l'édition de lien le symbole utilisé par la
    librairie hdf5 ne sera plus présent.
