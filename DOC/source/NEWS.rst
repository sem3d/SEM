.. -*- coding: utf-8; mode:rst -*-

Version 2020.03.0
-----------------

- Nombre de GLL enlev� des fichiers mat.dat et mater.in

- ngll defini dans input.spec (un seul valeur)

Version 2016.11.0
-----------------

- Loi materiau non-linéaire

- Les PML peuvent être définies sur plusieurs couches.

- Mise en données :

  - On peut maintenant spécifier des grilles cartésiennes de
    propriétés matériaux qui seront interpolées sur le maillage
    SEM. (uniformisation des champs aléatoires, modèles EarthChunk,
    PREM, etc...)

  - On peut définir des Ondes Planes, conditions de Neumann et de Dirichlet sur les surfaces
    du maillage

- Sorties :

  - (snapshots et capteurs) : ajout de l'energie d'onde P et S, et des
    tenseurs de contraintes. Possibilité de sortir les dérivées spatiales du
    déplacement pour les capteurs.

  - Sélection des types de sorties depuis le fichier ``input.spec``

  - On distingue des champs aux noeuds et aux mailles

- Optimisation du code :

  - restructuration du stockage interne des données pour permettre la vectorisation
    sur les architectures modernes

  - Compilation avec gfortran 6.2 / Intel 17

- Partitionneur :

  - optimisations, réécriture.

  - équilibrage de charge : prise en compte de coût différent selon le type de maille
    (solide, fluide, solide PML, fluide PML)

  - le mailleur automatique (cube) peu définir des PMLs sur plusieurs couches

  - lecture de maillage au format Abaqus, UNV, HDF5, GMSH

- Génération de champs aléatoire :

  - la génération est maintenant indépendante du lancement de SEM.

  - La mise en données s'appuie sur la lecture de champ cartésiens avec interpolation

  - La génération des champs est parallélisée

- SEM2D : partitionneur, introduction d'éléments de type DG, PML convolutionnelles.

- Documentation : mise à jour et compléments

Travaux en cours:

- Optimisation : parallélisation hybride OpenMP

- Compilation en mode simple précision

- PML Convolutionnelles (3D)

Statistiques :

539 fichiers modifiés,
104879 insertions de lignes,
26981 suppressions,
10 contributeurs.

Version 2015.02.1
-----------------

- Correction de bugs

Version 2015.02.0
-----------------

Ci-dessous sont recensées les nouveautés et modifications apportées
depuis la dernière version :

- Développement d'un partitionneur pour SEM 2D / C++ (LA)

- Entrées/sorties HDF5 pour SEM2D, sauvegarde de l'accélération dans les snapshots (LA)

- Nouvelle option (``dim``) obligatoire, permet :

  - la direction des pulses est maintenant un vecteur quelconque.

  - les vecteurs sont spécifiés en 2D pour Sem2D (``dim=2``) et non plus en X (Y=0) Z.

  - tenseur moment quelconque en 2D

- Nouveau format pour la spécification des capteurs (intégré à ``input.spec``)

- Model global EarthChunk (ML)

- Calcul de l'énergie totale pour Sem 2D (ST) (option ``output_total_energy``)

- Mailleur et SEM3D :

  - condition de type Dirichlet (P=0) pour les fluides (LG)

- Génération de champs aléatoires de propriétés pour les paramètres de Lamé (RC, ECP)

- Corrections de bug :

  - Compilation en mode MPI=off.

  - Position des récepteurs, récepteurs et sources dans les éléments à 27 noeuds.

  - Fin de calcul pour différents cas (Solide-Fluide, Solide-Fluide avec atténuation).

  - Calcul des coefficients utilisés pour l'atténuation. (LG)

  - Protection/reprise pour les cas Solide-Fluide.

  - Les paramètres d'intégration de Newmark sont toujours fixés, mais
    ne l'étaient pas pour les PMLs, on pouvait donc avoir un jeu de
    paramètres différents pour le milieu solide et le milieu
    solide+PML.

- Performance :

  - Regroupement des sorties par groupes de processeurs (option ``group_outputs``).

  - SEM2D, sauvegarde par blocs des capteurs.

  - Format HDF5 pour les capteurs.

  - Optimisation mailleur 3D pour les gros cas (LG).

  - Optimisation routine de calcul de l'atténuation.

Version 2013.12
---------------

- Stabilisation des éléments fluides et couplage fluide/solide.

- Nouvelles fonctions d'évolution temporelle des sources (square, tanh, ...)

- La limite maximum du nombre de processeurs gérés par le mailleur est
  portée à 8192 (1024 précédemment). La consommation mémoire excessive
  (en nombre de processeurs * nombre de mailles) est résolue. On a
  maintenant une consommation proportionnelle à : nombre de mailles x
  nombre de processeurs voisins.

- Un paramètre ``amplitude`` global pour toutes les fonctions temporelles est ajouté.

- La dépendance sur la librairie HDF5 est maintenant obligatoire.

- Le mailleur accepte un format HDF5 (d'un structure semblable au format UNV) en entrée.
  Cela permet de gérer de gros maillages, beaucoup trop longs à lire dans un format texte.

- Limitation du nombre de sorties texte du code pour passer des codes sur un grand nombre
  de processeurs.

- Bugfix: le paramètre ``mpml`` était ignoré.

- Bugfix: le mailleur ne libérait pas immédiatement les ressources
  HDF5, ce qui induisait des temps très longs de flush en fin de
  job pour les gros maillages.

  Ce temps étant décompté après l'épilogue MPI (ie après l'appel à ``MPI_Finalize``),
  le processus était considéré comme bloqué et tué avant la fin par ``mpirun``.

- Bugfix: les paramètres d'atténuation n'étaient pas correctement
  sauvegardés dans les fichiers de protection/reprise.

- Les sorties ont été mutualisées par groupe de processeurs, permettant d'avoir
  une sortie par noeud de calcul, plutôt qu'une sortie par processus MPI.

- Bugfix : pour les fluides, les excitations sont pondérées par :math:`\lambda` et plus :math:`\rho`.

- Fluide : on assure la continuité de :math:`\rho{}\phi` et non plus :math:`\phi` seul.

- Mailleur: on peut mailler un milieu stratifié simple.

- SEM2D : mutualisation de code, utilisation du nouveau format de fichier d'entrée commun avec SEM3D.

- Nouveaux champs en sortie des snapshots : pression et accélération.

- Optimisation : concerne le calcul des forces solides sans Acoef.

  Le calcul des dérivées spatiales a maintenant des cas particuliers
  pour ngll=5 et 7. Au delà de 10, l'ancienne méthode avec DGEMM
  optimisée (MKL) devient plus intéressante.

- Correction d'une fuite mémoire (à l'initialisation) dans l'allocation des capteurs.

Version 2013.04
---------------

Cette version résulte de l'intégration dans :program:`RegSEM.U` de :

- des modifications apportées par la version interne CEA,

- des éléments fluides développés dans une autre version issue de :program:`RegSEM.U`,

- de nouveaux développements destinés à simplifier l'utilisation et la
  maintenance du code.


Les nouveautés
~~~~~~~~~~~~~~

On liste ici les nouvelles fonctionnalités par rapport à :program:`RegSEM.U`.

Fonctionnalités du code :

- (:program:`SEM3D`) Introduction d'éléments de type fluide, avec couplage fluide-solide.

- Introduction d'un mécanisme d'amortissement sismique. On spécifie :math:`Q_\kappa`
  et :math:`Q_\mu` dans le fichier matériau. La bande de fréquence et le
  paramétrage du filtre sont déterminés par le fichier de configuration.

- Nouvelles formes d'onde pour les sources (Benchmark E2VP, Benchmark
  SPICE, sinus).

- Une variante des PML (MPML) avec son paramètre associé a été
  introduite. Ceci afin de régler des problèmes d'instabilités
  constatés sur certains cas.

- Un mode couplage optionnel avec un code externe.

- On peut maintenant faire des sorties snapshots partielles. Le fichier
  ``input.spec`` permet de décrire simplement une sélection de mailles
  à inclure dans les sorties.

Entrées/sorties :

- (MESH) Lecture des maillages au format unv.

- (:program:`SEM3D`, :program:`SEM2D`) Un nouveau format de fichier d'entrée (input.spec) :

  L'ancien format était très confus : une liste de valeurs lues de
  manière aveugle par les codes. Chaque code lisait ses paramètres
  dans un ordre pré-établi. Il était impossible de réutiliser un
  fichier de config d'une version à l'autre.

  Désormais les paramètres sont identifiés par des mots-clefs. Ainsi
  un paramètre inconnu est soit ignoré soit génère une erreur.

  Les sources sont décrites dans ce format.

- Les snapshots sont au format :program:`HDF5` :

  Le code génère en plus des fichiers :program:`HDF5`, un fichier XML (format
  XDMF) qui permet d'ouvrir directement les sorties dans :program:`Paraview` ou
  :program:`Ensight` (v10).

- Les maillages en entrée sont également au format :program:`HDF5` :

  Des problèmes de numérotation apparaissaient avec des gros maillages
  (utilisation du format ``I6`` pour les entiers). De plus, chacune des
  versions utilisait une variante subtile du même format texte (une
  ligne d'espacement pour l'un, un champ supplémentaire pour une
  autre...).

  Les identifiants sont maintenant des entiers 32 bits permettant de
  décrire 2 milliards de noeuds uniques, et le format utilise par
  défaut la compression gzip.

- Nouveau format pour le fichier des capteurs/traces :

  On a conservé le format de la version CEA, plus général. Dans une
  prochaine version ce fichier migrera vers un format semblable à
  celui de ``input.spec``.

- Le format des backups est désormais :program:`HDF5` (protection/reprise).

  Ce développement a été effectué pour faire passer un cas HPC. Le
  temps de création d'un backup pour ce cas est passé de 2H à 5min.

Optimisations :

- Optimisation des communications :

  L'algorithme d'échange inter-processeur a été entièrement revu pour
  utiliser des communications asynchrones. Il n'y a plus de risque
  d'interblocage occasionnel et les performances sont accrues.

- Optimisation de la consommation mémoire :

  Les mailles non-PML consommaient inutilement de la mémoire en
  stockant des pointeurs (non-alloués) vers des tableaux concernant
  uniquement les mailles PML.

  Une structure spécifique PML a été introduite. Celle-ci n'est
  allouée qu'au besoin uniquement pour les éléments contenant des PML.
  La mémoire utilisée est réduite à l'espace d'un seul pointeur par
  élément au lieu d'une dizaine.

- L'utilisation de la librairie :program:`HDF5` permet d'optimiser grandement les
  Entrées/Sorties pour les gros cas de calcul.


Autres :

- Améliorations du mailleur intégré :

  On utilise :program:`Metis` 5.x comme partitionneur. Ceci permet d'utiliser une
  topologie connectant toutes les mailles adjacentes (ayant au moins
  un vertex commun) contrairement à la version précédente qui ne
  considérait que les faces.

  Le mailleur génère ses maillages au format :program:`HDF5` attendu par SEM.

  De nombreuses optimisations et restructurations du code ont été
  effectuées accélérant le traitement.

- Introduction d'un répertoire de cas tests de non-régression et de
  benchmarks.

  Des cas d'exemples d'utilisation de :program:`SEM3D` se trouvent dans ``SEM3D/TESTS``.

- Compilation des sources avec :program:`CMake` :

  :program:`CMake` est un outil (comme autotools) permettant de générer des Makefiles.
  (voir :ref:`installation` ).

- Correction des FPML.

- (:program:`SEM3D`) : le code a été factorisé (suppression des duplications,
  réorganisations, simplifications) en plusieurs endroits.

