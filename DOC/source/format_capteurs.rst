.. -*- coding: utf-8; mode:rst -*-

Format du fichier de description des capteurs
=============================================

Exemple ::
  
  format = txt; # ou hdf5

  # 4 Capteurs specifies par leurs coordonnees
  # Les sorties sont nommees pts0001 .. pts0004
  capteur "pts" {
    type = points;
    count = 4;
    coords = "fichier_coords.txt";
    periode = 1;
  };
  # 50 Capteurs, sur une ligne specifiee par ses deux extremites.
  # Les sorties sont nommees line0001 .. line0050
  capteur "line" {
    type = line;
    count = 50;
    point0 = 0. 0. 0.;
    point1 = 0. 5000. 0.;
    periode = 1;
  };

  # Non implemente :
  capteur nom_capteur {
    coords = 0. 0. 0.;
    ligne = from 0. 0. 0. to 1. 1. 1. with 10; # Defini 10 capteurs entre les deux coords
    periode = 10; # default 1, intervalle entre 2 sauvegardes
    rayon = 10; # default 0 : pt de gauss le plus proche
    operation = avg; # kw: min, max, avg, gauss, interp
    grandeurs = vel acc pression; # possibles: vel acc pression depla deformation contrainte
  };
  ligne = from 0. 0. 0. to 1. 1. 1. with 10; # Defini 10 capteurs entre les deux coords
    periode = 10; # default 1, intervalle entre 2 sauvegardes
    rayon = 10; # default 0 : pt de gauss le plus proche
    operation = avg; # kw: min, max, avg, gauss, interp
    grandeurs = vel acc pression; # possibles: vel acc pression depla deformation contrainte
  };



Format du fichier de description des materiaux
==============================================

Exemple ::

  material N {
    # Version SEM
    domaine = fluide|solide;
    ngll = 5;

    vit 0 {  # 0 pour constante, 0 et 1 pour gradients
      pspeed = 2500;
      sspeed = 1500;
      Qp = 100;
      Qs = 100;
    };

    gradient {
      plane = ax+by+cz+d; # defini le plan 0, le plan 1 est implicitement defini parallele avec (d+1) comme offset
      z0 = 1000;
      z1 = 0;
    };

    elas 0 {
       # Definition par E,nu
    };

    compaction 0 {
    }

    hugoniot 0 {
    };

    plasticite 0 {
    };

    ngll=5;


  };
