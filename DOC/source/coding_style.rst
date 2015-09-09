.. -*- coding: utf-8; mode: rst -*-

=======================
Règles de programmation
=======================

Chaque développeur est invité, s'il le peut, à suivre les règles suivantes :

#. Au début de chaque fichier source, mettre (et compléter) cette en-tête : ::

    !>
    !! \source.f90
    !! \brief une description brève du contenu du fichier source
    !!
    !<

#. A la fin de chaque fichier source, mettre ce pied de page (afin que l'éditeur de texte sache indenter le code de façon homogène) : ::

    !! Local Variables:
    !! mode: f90
    !! show-trailing-whitespace: t
    !! f90-do-indent: 4
    !! f90-if-indent: 4
    !! f90-program-indent: 4
    !! f90-continuation-indent: 4
    !! End:
    !! vim: set sw=4 ts=8 et tw=80 smartindent : !!

   Si l'éditeur de texte ne reconnaît pas ces balises (eclipse), paramétrer l'éditeur pour qu'il remplace les tabulations par 4 espaces

#. Afin de limiter les fuites mémoires, lorsque l'on a besoin d'introduire un nouveau pointeur :

   #. Ajouter d'abord la déallocation du pointeur (typiquement dans deallocate_domain:deallocate_domain.f90)

   #. Ajouter ensuite l'allocation du pointeur (typiquement dans allocate_domain:allocate_domain.f90)

   #. Utiliser le pointeur (dans le code exécuté entre allocate_domain et deallocate_domain)

#. Vérifier qu'un pointeur est valide avant de l'utiliser. Si le pointeur est invalide, arrêter le programme et afficher un message d'erreur. Par exemple, en fortran : ::

     if ( .not. allocated ( ptr ) ) stop 'ERREUR - NomDeLaSubroutine : données invalides'

#. En fortran, utiliser des fonctions / routines préférentiellement placées dans des modules. Le compilateur fortran ne fait pas de vérification sur les arguments (type, nombre) lorsqu'il compile une fonction / routine qui n'est pas placée dans un module

#. En fortran, utiliser ``implicit none`` dans les fonctions / routines

#. Ne pas hésiter à couper un module fortran en 2 (ou plus) parties s'il devient trop grand

#. Choisir des noms clairs et explicites pour nommer les fonctions / routines

#. Limiter la taille (2 écrans maximum) de chaque fonction / routine. Re-découper une fonction / routine si elle devient trop grande

#. Initialiser chaque variable (et pointeur) avant de l'utiliser

#. Factoriser (ré-utiliser) le code le plus possible quand cela est possible

#. Après chaque développement, ajouter un test à la suite de tests pour valider le développement

#. Après chaque développement, rejouer la suite de tests fournie avec le code et vérifier que les tests passent
