.. -*- coding: utf-8; mode:rst -*-

Format du fichier de description des capteurs
=============================================

Exemple ::
  
  format = txt; # ou hdf5

  # 4 Capteurs specifiés par leurs coordonnées
  # Les sorties sont nommées pts0001 .. pts0004
  capteur "pts" {
    type = points;
    coords = "fichier_coords.txt";
    period = 1;
  };
  # 50 Capteurs, sur une ligne spécifiée par ses deux extrémités.
  # Les sorties sont nommées line0001 .. line0050
  capteur "line" {
    type = line;
    counti = 50;
    point0 = 0. 0. 0.;
    point1 = 0. 5000. 0.;
    period = 1;
  };


Format du fichier de description des matériaux
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
      plane = ax+by+cz+d; # défini le plan 0, le plan 1 est implicitement défini parallèle avec (d+1) comme offset
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
