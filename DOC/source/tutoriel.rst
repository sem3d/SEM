.. -*- mode:rst; coding: utf-8 -*-

======================
Prise en main de SEM3D
======================

Introduction
============

Ce chapitre présente la chaine de calcul SEM, insistant plus
particulièrement sur SEM3D.

La chaîne logicielle SEM contient le code de simulation (``sem2d.exe``
et ``sem3d.exe``, ainsi qu'un outil de partitionnement et préparation
de maillage ``mesher``).

L'outil ``mesher`` permet de partitionner des maillages complexes au
format abaqus, UNV ou encore un format spécifique simple (équivalent
du format UNV dans un fichier HDF5).

Outre la paramétrisation du code SEM3D, la plus grosse difficulté dans
l'utilisation du code concerne la création des maillages, qui doit être
effectuée par des outils spécialisés comme Cubit. Pour des maillages
simples : stratifiés, avec topographie, on peut utiliser la suite
d'outils *meshtools*.

Les équations du mouvement
--------------------------

SEM résout la propagation d'onde élastique dans un milieu décrit par l'équation
locale du mouvement reliant le déplacement :math:`u` en chaque point matériel, les
contraintes :math:`\sigma` et les forces extérieures :math:`\vec{f}` :

.. math::

   \rho \frac{\partial^2 u}{\partial t^2} = \nabla.\sigma + \vec{f}

Avec en élasticité linéaire : :math:`\sigma=C:\nabla{}u`, où :math:`C` est le
tenseur élastique d'ordre 4.

Pour l'instant les milieux de propagations décrits dans SEM sont
considérés isotropes.  Le code est prévu pour gérer les milieux
anisotropes, mais il n'existe pas de manière simple de gérer la mise
en données.

Formulation éléments finis
--------------------------

SEM est un code éléments finis, basé sur une formulation spectrale
(d'où son nom). Le champ de déplacement :math:`u` est décrit dans
chaque élément, ou maille, sur une base de polynômes de Lagrange
d'ordre N (N défini comme paramètre).

Pour obtenir une convergence spectrale, ces polynômes de Lagrange sont
définis sur les points de Gauss-Lobato-Legendre (GLL) de chaque
éléments (voir :ref:`fig-gll`).

.. _fig-gll:

.. figure:: ../figures/elem_gll_2d.png
   :scale: 40%
   :align: center

   Position des points de Gauss-Lobato-Legendre sur un élément 2D d'ordre 9
   

Les éléments traités sont des quadrangles en 2D et des Hexaèdres en
3D. Si :math:`\Phi_i(x)` est le polynôme de Lagrange valant 1 au point
GLL :math:`x_i`, le champ de déplacement :math:`u(x,y,z)` de l'élément
s'exprime sur la base tensorisée :math:`\Phi \otimes \Phi \otimes
\Phi` :

.. math::

   u(x,y,z) = \sum_{i,j,k} U_{i,j,k}.\Phi_i(x).\Phi_j(y).\Phi_k(z)

Ainsi, sur un élément d'ordre 5, la composante *x* du champ de
déplacement, est décrite par un vecteur de 125 éléments
:math:`U_{i,j,k}` .

La figure :ref:`pollag` montre la forme des polynôme de Lagrange d'ordre 9, la base tensorisée
de dimension 2D est représentée :ref:`fig-ref-2d`

.. _pollag:

.. figure:: ../figures/lagrange_9.png
   :scale: 40%
   :align: center

   Polynômes de Lagrange d'ordre 9.


.. _fig-ref-2d:

.. figure:: ../figures/shape_fct_9x9.png
   :scale: 60%
   :align: center

   Quelques fonctions de forme d'un élément 2D d'ordre 9

Le champ de déplacement est continu entre deux éléments adjacents, et
le maillage géré par SEM doit être conforme (toutes mailles se
touchant ont en commun un sommet, ou une arête ou une face
complète). De plus l'ordre en X, Y ou Z de chaque maille doit assurer
la conformité au niveau des points GLL en commun.

Enfin, dans ce qui précède, on a présenté la formulation de manière
simplifiée, sur des mailles cubiques, alignées en *x,y,z*. En pratique
SEM peut gérer un maillage hexaédrique quelconque (mais conforme),
composé de mailles parallélépipédiques. Dans chaque élément le code se
ramène à une base locale :math:`\phi_i(x).\phi_j(y).\phi_k(z)` par
changement de variable de la fonction de base :math:`\Phi` depuis la
maille vers un élément de référence sur le segment [-1.,1.] .

Enfin une des originalités de la méthode, provient du choix de
quadrature pour l'évaluation numérique des intégrales apparaissant
dans la formulation élément finis.

On passe d'abord à la formulation faible en multipliant l'équation
locale par une fonction quelconque *w* et en intégrant (produit
scalaire dans :math:`\mathcal{L}^2`)

.. math::

   \int w.\rho \frac{\partial^2 u}{\partial t^2}\vec{dx} = \int w.(\nabla.(C:\nabla{}u) + \vec{f}).\vec{dx}

En exprimant *w* et *u* sur la même base discrète
:math:`\Phi_i(x,y,z)` (ici *i* indexe **toutes** les fonctions de base
de tous les éléments).

.. math::

   \sum_{i,j} w_i.\rho \frac{\partial^2 u_j}{\partial t^2}\int \Phi_i\Phi_j \vec{dx} = 
     \sum w_i.u_j.(\nabla.(C:\nabla{}\Phi_j) + f_j\Phi_j).\vec{dx}

Cette dernière équation apparaît alors sous la forme classique de
l'approximation de Galerkin : :math:`a(u,w) = f(w)` avec :math:`a` une
forme bilinéaire.

Sans aller jusqu'au bout des développements, on voit qu'il apparaît une
matrice :math:`M_{i,j}=\int \Phi_i\Phi_j\vec{dx}`, que l'on doit
inverser si on veut obtenir une expression de :math:`\frac{\partial^2
u_j}{\partial t^2}` .

Les produits scalaires entre fonctions :math:`\Phi_i` qui ne partagent
pas le même élément support sont nuls par construction. Mais au sein
d'un éléments, les polynômes de Lagrange ne sont pas orthogonaux. La
méthode SEM utilise astucieusement une quadrature basée sur les mêmes
points de Gauss que les noeuds de définitions des fonctions de
base. Cela introduit bien sûr une approximation de l'intégrale, mais
le résultat est que le produit scalaire discret utilisé rend
orthogonale les fonctions :math:`\Phi_i` ayant le même élément
support.


Conditions de bord
------------------

La condition naturelle d'un bord en élément fini est d'être une
surface libre, donc réfléchissante pour les ondes. Pour simuler des
milieux ouverts, SEM implémente un type d'élément dit PML (Perfectly
Matched Layer) pour simuler un milieu ouvert infini en bordure d'un
domaine.

Intégration temporelle
----------------------

Le schéma d'intégration est un schéma de Newmark explicite.


Le pas de temps d'intégration dans SEM est calculé automatiquement à
partir du nombre de Courant :math:`\mathcal{C}<1` (paramètre de configuration) selon :

.. math::

   \Delta t = \mathcal{C} \frac{\min \Delta{x}}{\max Velocity}

Attention:

   Des mailles trop petites, ou des vitesses de propagation trop
   importantes vont faire chuter le pas de temps.

Résolution spatiale
-------------------

Le maillage doit également être suffisement résolu pour capturer les
fréquences spatiales du signal que l'on veut propager. On considère
que 10 points GLL par longueur d'onde sont suffisant.

Augmenter l'ordre des éléments est donc un moyen d'obtenir une
résolution spatiale correcte avec un maillage donné. La convergence
spatiale étant rapide, augmenter l'ordre devrait permettre de baisser
le nombre de points par longueur d'onde nécessaire, mais cela augmente
doublement les coûts de calcul :

- la complexité est en :math:`N^3` par points GLL,

- le pas de temps est proportionnel à :math:`\frac{1}{\min \Delta x}`,
  le pas d'espace :math:`\min \Delta x` diminuant avec l'ordre des
  éléments (On voit sur :ref:`fig-gll` comment les points de Gauss se
  ressèrent vers les bords avec l'augmentation de l'ordre.

Atténuation
-----------

Un mécanisme d'atténuation sismique des ondes P et S est implémenté,
sous forme d'une série de filtres répartis sur une bande de
fréquence. (voir [KOM98]_)


Description des sorties
-----------------------

Les résultats de simulation peuvent être obtenus sous deux formes :

- Des instantanés (*snapshot*) des champs obtenus sur tous les points GLL, ou sur
  un sous-partie, à une fréquence données. Ces sorties sont en général
  assez lourdes et ne peuvent être trop fréquentes.

- Des sorties *capteurs*, pour un ou plusieurs points du maillage, on
  sort les valeurs du champ toutes les N itérations de calcul.

Les champs disponibles sont :

============= ====== ======== ========
Champ         Milieu Snapshot Capteurs
============= ====== ======== ========
Déplacement   S      Oui      Oui
Vitesse       S/F    Oui      Oui
Accélération  S/F    Oui      Non
Pression      S/F    Oui      Non
============= ====== ======== ========


Pour les instantanés, il existe un mécanisme de sélection de mailles
qui permet de ne sauvegarder qu'une partie du maillage. Cependant on
ne peut sélectionner que des mailles complètes (donc avec tous ses
points GLL), et pour l'instant, on ne peut pas, sauf en
post-traitement, réinterpoler les fonctions de formes sur un maillage
plus grossier.

Présentation des outils
=======================

Deux exécutables sont impliqués directement dans l'utilisation de SEM :

- ``mesher`` et ``sem3d.exe`` pour le cas 3D,

- ``sem2d.exe`` pour le cas 2D, il n'existe pas encore d'outil de
  partitionnement simple à utiliser.

``mesher`` transforme un maillage d'entrée en un maillage partitionné
utilisable par SEM. On peut lui fournir différents formats :

- Un maillage au format *Abacus* (d'extension ``.aba``)

- Un maillage au format *UNV*, (aussi connu sous le nom *IDEAS*)
  d'extension ``.unv``, contenant des hexaèdre pour la 3D.

- Un maillage au format *HDF5*, spécifique, dont la structure est
  décrite en détail dans _`Format HDF5`, contenant 3 tables :

  - ``Elements`` : un tableau de NE x 8 entiers de 0 à (NN-1) faisant
    référence aux noeuds.

  - ``Nodes`` ; un tableau de NN x 3 de rééls, les coordonnées des
    noeuds

  - ``Mat`` : un tableau de NE entiers, contenant le numéro matériau à
    associer à chaque maille.

- Le quatrième format est simplement la description d'un maillage
  cartésien, pour lequel on entre manuellement les coordonnées et la
  subdivision de la grille souhaitée.


L'outil mailleur, en plus de ses entrées en ligne de commande,
s'appuie sur un fichier externe ``mat.dat``, donnant quelques
informations sur le maillage à générer : nombre de matériaus, présence
d'éléments PML, type de matériau (solide ou fluide).


Préparation d'un cas de calcul
------------------------------

Pour lancer un calcul SEM, il faut se placer dans le répertoire du cas et y placer
les fichiers nécéssaires à son exécution. L'arborescence doit être la suivante ::

  CAS/
  |- input.spec
  |- material.input
  |- sem/
  |  |- mesh4spec.0000
  |  |- ...
  |  |- mesh4spec.NNNN
  |- capteurs.dat

input.spec:

  Ce fichier contient la configuration du code :
  - paramètres d'intégration temporelle, temps physique du calcul,
  - description de la ou des sources,
  - description des sorties capteurs,
  - description des sorties snapshots.

material.input:

  Ce fichier contient la description de chaque matériau : :math:`\rho, V_p, V_s`, un nombre
  de points GLL par direction de la maille de référence.

capteurs.dat

  Contient une description des sorties capteurs souhaitées.

Le fichier ``input.spec`` est décrit en détail dans la section
_`Description des paramètres de SEM3D`.

Des exemples de fichiers `material.input` et `capteurs.dat` sont
disponibles dans les tests du code. Ces derniers sont de simples
tables de paramètres.


Exemples de modélisation avec SEM3D
===================================


Maillage uniforme avec PML
--------------------------

On commence par un premier exemple de grille cartésienne avec une
source ponctuelle.

Le fichier ``mat.dat`` doit contenir (les commentaires, après le *#*
sont facultatifs) ::

  1  # number of non PML materials
  F  # Milieu stratifié F: non T: oui
  1  # PMLs? 0: no, 1: yes
  1 1  # PMLs on top? at the bottom? (0: no, 1: yes)
  S

On lance l'exécutable ``mesher``, et on lui indique les informations
suivantes :

- Nombre de procésseurs : 4

- Construction du modèle matériaux et maillage : 1 (Oui)

- Choix d'une grille : 1 (On the fly : sur la mouche)

- Saisie des coordonnées et taille de maille : 

  - X : -100, 500
  - Y : -100, 500
  - Z : -100, 500

  - DX, DY, DZ : 50

- Choix de 8 noeuds par maille : 1 (Les mailles quadratiques à 27
  noeuds sont en développement)

L'outil va alors générer 4 fichiers nommés ``mesh4spec.000N.h5``
(N=0,1,2,3) contenant les maillages et informations de communication
des 4 partitions.

Lancement du cas
----------------

Il faut d'abord préparer le répertoire du CAS : y copier les fichiers
``input.spec``, ``material.input``, ``capteurs.dat``, et placer les fichiers
``mesh4spec.NNNN`` dans le sous-répertoire ``sem/``.


