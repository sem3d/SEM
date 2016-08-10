.. -*- coding: utf-8 -*-

=============================================
Annexe: Description des Datasets msh de Gmesh
=============================================


source: http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-1_002e0

#. Générer un format `msh` dans Gmesh 
   
   Les formats de fichier `msh` peut être lu par les mailleurs disponible dans SEM. 
   Cependant, le format de fichier compatible est celui celui sauvegardé en `ASCII`
   en laissant également toutes les options non `cochée` comme illustré sur la figure
   ci-dessous :
      
       .. _fig-mesh-msh:

       .. figure:: ../figures/save_gmesh.png   
          :scale: 50
          :align: center
  
   Une deuxième possibité consister à sauvegader via le menu-déroulante dans `Gmesh` 
   (cf figure ci-dessous).
    
       .. _fig-mesh-msh:
         
       .. figure:: ../figures/save_msh.png
          :scale: 50
          :align: center
   
   Ces deux options de sauvegarde permettent de sortire dans le fichier `msh` des éléments
   finis de type hexaèdrique et/ou quadratique, qui sont les seuls types d'éléments finis 
   manipulable dans SEM.       
   
   Pour une simulation 3D, les éléments finis de type quadratique servent à définir
   les différentes surfaces qui sont évoquée dans la section surface.
        
#. Description du fichier de sortie
      
   Chaque bloc d'information (ensemble de données) est délimité par `\$Nom_bloc` pour indiquer le
   début et `\$EndNom_block` pour marquer la fin du bloc. Le lecteur de SEM ignorera toute ligne
   commançant par `\$`.
   
   
   **Dataset \$MeshFormat**
   
   Il s'agit du premier bloc d'information apparaissant dans le fichier `.msh` livrant les information générales
   sur ce dernier. Ce bloc est rganisé de la façon suivante : ::
   
     $MeshFormat
      version-num file-type data-size
     $EndMeshFormat
   
   `version-num` est un réel indiquant la version de Gmesh utilisé
   
   `file-type` vaut 0 pour un fichier texte.
   
   `data-size` indique le nombre de décimales significatives
   
   Ci-dessous un exmaple de ce bloc ::
    
     $MeshFormat
      2.2 0 8
     $EndMeshFormat
     
   
   **Dataset \$Nodes** ::
     
   Il s'agit du deuxième bloc par ordre d'apparition dans le fichier donnant des information sur les noeuds du maillage. Il fournit 
   comme informations essentielles le nombre de noeuds et la liste des coordonnées de ces derniers. ::
     
     $Nodes
      number-of-nodes
      node-num x y z
      ...
     $EndNodes  
   
   `number-of-nodes`: indique le nombre de noeuds générés pour faire le maillage.
   
   `node-num`       : numérotation global des noeuds
   
   `x y z`          : coordonnées globales x, y, z de chaque noeud
   
   Ci-dessous un example ::
   
     $Nodes
      7
      1 -1 -1 0
      2 1 -1 0
      3 1 1 0
      4 -1 1 0
      5 -1 -1 1
      6 1 -1 1
      7 1 1 1
     $EndNodes
   
   **Dataset \$Elements** ::
     
   Il s'agit du troisième bloc dans le fichier `.msh` qui définit les différents types d'éléments utilser
   et les différents groupes phyiques auquel ils appartiennent. Ce bloc est organiser de la façon suivante : ::
   
     $Elements
       number-of-element
       elm-num elmt-type number-of-tags list-tags list-of-nodes
       ....
   
     $EndElements
   
   `number-of-element`: indique le nombre d'éléments présents dans présent dans le maillage
   
   `elm-num`          : le numéros global attribué à un élément fini
   
   `elmt-type`        : indique le type de l'élément finis (géormétrie de la maille) `elm-num`

                       - `elmt-type` = 3 pour les éléments finis quadrangle à 4 noeuds

                       - `elmt-type` = 5 pour les éléments finis hexaédrique à 6 faces
   
   `number-of-tags`   : indique le nombre de paramètres associés à l'élément
   
   `list-tags`        : indiques les paramères associés à l'éléments et est de taille `number-of-tags`
   
   `list-of-nodes`    : donne la liste ordonnée des références des noeuds constituant l'élément.
   
   Ci-dessous un example de ce bloc : ::
   
      $Elements
      343
      1 5 2 1 1 24 40 5 1 41 57 25 9
      2 5 2 1 1 23 39 40 24 42 58 57 41
      ....
   
      18 5 2 1 1 51 67 66 50 55 71 70 54
      ...
      ...
      $EndElements
   
   **Example de maillage au format msh** ::
   
       $MeshFormat
       2.2 0 8
       $EndMeshFormat
       $Nodes
       75
       1 -1 -1 0
       2 1 -1 0
       3 1 1 0
       4 -1 1 0
       5 -1 -1 1
       6 1 -1 1
       7 1 1 1
       8 -1 1 1
       9 -0.5000000000011573 -1 0
       10 -2.752797989558076e-12 -1 0
       11 0.4999999999986235 -1 0
       12 1 -0.5000000000011573 0
       13 1 -2.752797989558076e-12 0
       14 1 0.4999999999986235 0
       15 0.5000000000011573 1 0
       16 2.752797989558076e-12 1 0
       17 -0.4999999999986235 1 0
       18 -1 0.5000000000011573 0
       19 -1 2.752797989558076e-12 0
       20 -1 -0.4999999999986235 0
       21 -0.5000000000011573 -1 1
       22 -2.752797989558076e-12 -1 1
       23 0.4999999999986235 -1 1
       24 1 -0.5000000000011573 1
       25 1 -2.752797989558076e-12 1
       26 1 0.4999999999986235 1
       27 0.5000000000011573 1 1
       28 2.752797989558076e-12 1 1
       29 -0.4999999999986235 1 1
       30 -1 0.5000000000011573 1
       31 -1 2.752797989558076e-12 1
       32 -1 -0.4999999999986235 1
       33 -1 -1 0.4999999999986718
       34 1 -1 0.4999999999986718
       35 1 1 0.4999999999986718
       36 -1 1 0.4999999999986718
       37 -7.578315348306864e-24 0 0
       38 -1.376398994782827e-12 -0.5 0
       39 -0.5 1.376398994779038e-12 0
       40 -0.5000000000005786 -0.4999999999993118 0
       41 1.376398994775249e-12 0.5 0
       42 -0.4999999999993118 0.5000000000005786 0
       43 0.5 -1.376398994779038e-12 0
       44 0.4999999999993118 -0.5000000000005786 0
       45 0.5000000000005786 0.4999999999993118 0
       46 -2.752909011860538e-12 -1 0.5
       47 -0.5000000000011573 -1 0.4999999999993359
       48 0.4999999999986238 -1 0.4999999999993359
       49 1 -2.752909011860538e-12 0.5
       50 1 -0.5000000000011573 0.4999999999993359
       51 1 0.4999999999986238 0.4999999999993359
       52 2.752909011860538e-12 1 0.5
       53 0.5000000000011573 1 0.4999999999993359
       54 -0.4999999999986238 1 0.4999999999993359
       55 -1 2.752909011860538e-12 0.5
       56 -1 0.5000000000011573 0.4999999999993359
       57 -1 -0.4999999999986238 0.4999999999993359
       58 -7.578315348306864e-24 0 1
       59 -1.376398994782827e-12 -0.5 1
       60 -0.5 1.376398994779038e-12 1
       61 -0.5000000000005786 -0.4999999999993118 1
       62 1.376398994775249e-12 0.5 1
       63 -0.4999999999993118 0.5000000000005786 1
       64 0.5 -1.376398994779038e-12 1
       65 0.4999999999993118 -0.5000000000005786 1
       66 0.5000000000005786 0.4999999999993118 1
       67 -7.578315348306864e-24 0 0.5
       68 -1.376454505934058e-12 -0.5 0.5
       69 -0.5 1.376454505930269e-12 0.5
       70 -0.5000000000005786 -0.4999999999993117 0.4999999999996679
       71 1.37645450592648e-12 0.5 0.5
       72 -0.4999999999993118 0.5000000000005786 0.4999999999996679
       73 0.5 -1.376454505930269e-12 0.5
       74 0.4999999999993118 -0.5000000000005785 0.4999999999996679
       75 0.5000000000005785 0.4999999999993119 0.4999999999996679
       $EndNodes
       $Elements
       56
       1 3 2 28 25 4 18 56 36
       2 3 2 28 25 18 19 55 56
       3 3 2 28 25 56 55 31 30
       4 3 2 28 25 36 56 30 8
       5 3 2 28 25 19 20 57 55
       6 3 2 28 25 20 1 33 57
       7 3 2 28 25 57 33 5 32
       8 3 2 28 25 55 57 32 31
       9 3 2 27 26 5 21 61 32
       10 3 2 27 26 21 22 59 61
       11 3 2 27 26 61 59 58 60
       12 3 2 27 26 32 61 60 31
       13 3 2 27 26 31 60 63 30
       14 3 2 27 26 60 58 62 63
       15 3 2 27 26 63 62 28 29
       16 3 2 27 26 30 63 29 8
       17 3 2 27 26 22 23 65 59
       18 3 2 27 26 23 6 24 65
       19 3 2 27 26 65 24 25 64
       20 3 2 27 26 59 65 64 58
       21 3 2 27 26 58 64 66 62
       22 3 2 27 26 64 25 26 66
       23 3 2 27 26 66 26 7 27
       24 3 2 27 26 62 66 27 28
       25 5 2 1 1 1 9 40 20 33 47 70 57
       26 5 2 1 1 33 47 70 57 5 21 61 32
       27 5 2 1 1 9 10 38 40 47 46 68 70
       28 5 2 1 1 47 46 68 70 21 22 59 61
       29 5 2 1 1 20 40 39 19 57 70 69 55
       30 5 2 1 1 57 70 69 55 32 61 60 31
       31 5 2 1 1 40 38 37 39 70 68 67 69
       32 5 2 1 1 70 68 67 69 61 59 58 60
       33 5 2 1 1 19 39 42 18 55 69 72 56
       34 5 2 1 1 55 69 72 56 31 60 63 30
       35 5 2 1 1 39 37 41 42 69 67 71 72
       36 5 2 1 1 69 67 71 72 60 58 62 63
       37 5 2 1 1 18 42 17 4 56 72 54 36
       38 5 2 1 1 56 72 54 36 30 63 29 8
       39 5 2 1 1 42 41 16 17 72 71 52 54
       40 5 2 1 1 72 71 52 54 63 62 28 29
       41 5 2 1 1 10 11 44 38 46 48 74 68
       42 5 2 1 1 46 48 74 68 22 23 65 59
       43 5 2 1 1 11 2 12 44 48 34 50 74
       44 5 2 1 1 48 34 50 74 23 6 24 65
       45 5 2 1 1 38 44 43 37 68 74 73 67
       46 5 2 1 1 68 74 73 67 59 65 64 58
       47 5 2 1 1 44 12 13 43 74 50 49 73
       48 5 2 1 1 74 50 49 73 65 24 25 64
       49 5 2 1 1 37 43 45 41 67 73 75 71
       50 5 2 1 1 67 73 75 71 58 64 66 62
       51 5 2 1 1 43 13 14 45 73 49 51 75
       52 5 2 1 1 73 49 51 75 64 25 26 66
       53 5 2 1 1 41 45 15 16 71 75 53 52
       54 5 2 1 1 71 75 53 52 62 66 27 28
       55 5 2 1 1 45 14 3 15 75 51 35 53
       56 5 2 1 1 75 51 35 53 66 26 7 27
       $EndElements
          
      
        
