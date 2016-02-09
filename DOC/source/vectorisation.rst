.. -*- mode: rst; coding: utf-8 -*-

=================================
Note relatives à la vectorisation
=================================

Pour maximiser l'adéquation entre structures de données (membre d'un module) et architecture (vectorisation possible ou non), des macros sont définies dans le fichier index.h : il existe une macro (mc_) par structure (mb_).

Les règles à suivre pour utiliser ces macros sont :

- allocation de la strucutre : utiliser mc_ (mais pas mb_).

- accès (lecture, écriture) aux éléments de la structure : utiliser mc_ (mais pas mb_).

- désallocation : utiliser mb_ (mais pas mc_).
