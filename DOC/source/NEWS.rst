.. -*- coding: utf-8; mode:rst -*-

NEWS
====

Version 2013.02
---------------

Cette version résulte de l'intégration de plusieurs version de SEM3D,
RegSEM.U, SEM-CEA.

Les nouveautés par rapport à toutes les versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Un nouveau format de fichier d'entrée (input.spec) :

  L'ancien format était très confu : une liste de valeurs lues de
  manière aveugle par les codes. Chaque code lisait ses paramètres
  dans un ordre pré-établi. Il était impossible de réutiliser un
  fichier de config d'une version à l'autre.

  Désormais les paramètres sont identifiés par des mots-clefs. Ainsi
  un paramètre inconnu est soit ignoré soit génère une erreur.

  Les sources sont décrites dans ce format.

- Les snapshots sont au format HDF5 :

  Le code génère en plus des fichiers HDF5, un fichier XML (format
  XDMF) qui permet d'ouvrir directement les sorties dans Paraview ou
  Ensight (v10).

- Les maillages en entrée sont également au format HDF5 :

  Des problèmes de numérotation apparaissaient avec des gros maillage
  (utilisation format I6 pour les entiers). De plus, chacune des
  versions utilisait une variante subtile du même format texte (une
  ligne d'espacement pour l'un, un champ supplémentaire pour une
  autre...).

  Les identifiants sont maintenant des entiers 32 bits permettant de
  décrire 2 milliards de noeuds uniques, et le format utilise par
  défaut la compression gzip.

- Optimisation des communications :

  L'algorithme d'échange inter-processeur à été entièrement revu pour
  utiliser des communications asynchrones. Plus de risque
  d'interblocage et des performances accrues.

- Optimisation de la consomation mémoire :

  Les mailles non-PML consommaient inutilement de la mémoire en
  stockant des pointeurs (non-alloués) vers des tableaux concernant
  uniquement les mailles PML.

  Une structure spécifique PML a été introduite qui n'est allouée
  qu'au besoin, consommant ainsi l'espace d'un seul pointeur au lieu
  d'une dizaine.

- Introduction d'éléments de type fluide, avec couplage fluide solide.

- Corrections de bugs :

  - Il manquait une équation dans le calcul des PML classiques.

- Améliorations du mailleur intégré :

  On utilise metis 5.x comme partitionneur. Ceci permet d'utiliser une
  topologie connectant toutes les mailles adjacentes (ayant au moins
  un vertex commun) contrairement à la version précédente qui ne
  considérait que les faces.

  Le mailleur génère ses maillages au format HDF5 attendu par SEM.

  De nombreuses optimisations on été efféctuées accélérant le
  traitement.

- Introduction d'un répertoire de cas tests de non-regression et de
  benchmarks.

  Les tests SEM3D se trouvent dans ``SEM3D/TESTS``

- Compilation des sources avec CMake :

  CMake est un outil (comme autotools) permettant de générer des Makefiles.
  (voir `DOC/INSTALLATION.pdf`_ )

Les nouveautés de cette version par rapport à RegSEM.U
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Nouveau format pour le fichier des capteurs/traces :

  On a conservé le format de la version CEA, plus général. Dans une
  prochaine version ce fichier migrera vers un format semblable à
  celui de ``input.spec``

- Nouvelles formes d'onde pour les sources. (Benchmark E2VP, Benchmark
  SPICE, sinus)

- gradient

- separation sources, gradients, etc...

- Format des backups en HDF5 (protection/reprise)

  Les fichiers de backup (ou protection/reprise) sont également en HDF5.

  Ce développement à été effectué pour faire passer un cas HPC. Le
  temps de création d'un backup pour ce cas est passée de 2H à 5min.

- Correction du calcul des mécanisme d'amortissement.

- Un mode couplage optionnel avec un code externe (pour l'instant
  Mka3D).

- Une variante des PML (MPML) avec son paramètre associé a été
  introduite. Ceci afin de régler des problèmes d'instabilités
  constatés sur certains cas.



TODO
~~~~

XXX: A décrire

- modele de source

- description condition Neumann

- description des capteurs (mot-clefs)

- anisotropie


Notes importantes
~~~~~~~~~~~~~~~~~

- Le schéma en temps à été simplifié (Les paramètres beta/gamma de
  l'algorithme de Newmark ne sont plus modifiables.

  Ils pourront être réintroduit une fois réglé le problème de
  synchronisation avec les forces de couplage externes.

- Bien que les deux méthodes continuent de coéxister, le calcul des
  forces utilisant le tableau Acoeff a été désactivé dans cette
  version. Le code est plus lisible mais moins rapide.

  On étudiera comment obtenir le meilleur des deux mondes dans une
  prochaine version.

