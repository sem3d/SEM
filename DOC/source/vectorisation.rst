.. -*- mode: rst; coding: utf-8 -*-

================================
Note relative à la vectorisation
================================

Pour maximiser l'adéquation entre structures de données (membre d'un module) et architecture (vectorisation possible ou non), des macros sont définies dans le fichier index.h : il existe une macro (suffixée par ``_``) par structure (préfixée par ``m_``).

Les règles à suivre pour utiliser ces macros sont :

- allocation et accès (lecture, écriture) : utiliser ``*_`` (mais pas ``m_*``). Pour pouvoir utiliser les macros il faut respecter la casse de leurs noms (le préprocesseur qui les substitue est sensible à la casse).

- désallocation : utiliser ``m_*`` (mais pas ``*_``).

Dans la configuration CMake, il est préférable d'utiliser OPT_VEC=OFF (ou OPT_VEC=ON avec soit OPT_VEC_LOW_RAM=ON, soit OPT_VEC_MED_RAM=ON) sur un portable, et, d'utiliser OPT_VEC_LOW_RAM=OPT_VEC_MED_RAM=OFF sur des machines vectorielles (AVX512).
