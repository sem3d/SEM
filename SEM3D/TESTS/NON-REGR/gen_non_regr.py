# -*- mode: python; coding: utf-8 -*-
# Outil pour la génération de cas de non-régression

import sys
import os

template = sys.argv[1]
outputdir = sys.argv[2]

MATERIAL = [
    ("S", 6300., 2500., 2800., 5, 0., 0.), # S/F Vp Vs Rho ngll Qp Qs
    ("F", 1300.,    0., 1000., 5, 0., 0.),
]

MATDEF="%(t)s %(vp)s %(vs)s %(rho)s  %(ngll)s  %(ngll)s  %(ngll)s 0.1 %(qp)s %(qs)s"
TXTLAYER="Epaisseur %(laythick)s m %(ndiv)s mailles"
LAYER="%(laythick)s %(ndiv)s # Layer $(laynum)s"

subst = {
    "name" : "Test_name",
    "save_snap" : "true",
    "snap_interval" : 0.01,
    "prorep" : "false",
    "restart_iter" : 0,

# Nb procs
    "nprocs" : 1,

# Materials
    "nmat" : 1,
    "matdef" : "S  6300.00  2500.00   2800. 5 5 5  0.000005 0. 0.",
# Dimensions du domaine
    "xmin" : -100,
    "xmax" :  500,
    "ymin" : -100,
    "ymax" :  500,
    "zmin" : -100,
    "zmax" :  500,
    "xstep":   50,
    "ystep":   50,
# Description couches (mat.dat)
    "nblayers" : 1,
    "layers"  : "600 12 # Layer 1",
    "txtlayer" : "Epaisseur 600m 12 mailles",
    "pmla" : 1,
    "pmlt" : 1,
    "pmlb" : 1,
    "pmlgll" : 5,
    "ctlpt" : 1,
# Position des tranches de snapshot
    "x0"   :  100,
    "x1"   :  150,
    "y0"   :  100,
    "y1"   :  150,
    "z0"   :  100,
    "z1"   :  150,
# Position de la source
    "xs"   : 25.,
    "ys"   : 25.,
    "zs"   : 25.,
    "stype" : "impulse",
    "sfunc" : "ricker",
}


FILES=["input.spec",
       "mat.dat",
       "mater.in",
       "mesh.input",
       "README",
]


def gen_files(template, src, dest, subst):
    for fname in FILES:
        data = open(pjoin(src, fname), "r").read()
        data = data % subst
        output = open(pjoin(dest, fname), "w")
        output.write(data)


GENLAYERS = [
    ("S__p%03(nprocs)d_g%(nglls)s", [ (600, 12, MATERIAL[0]) ]),
    ("F__p%03(nprocs)d_g%(nglls)s", [ (600, 12, MATERIAL[1]) ]),
    ("SF_p%03(nprocs)d_g%(nglls)s", [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ]),
]

GENPROCS = [ 1, 2, 4 ]

GENGLL = [ 5, 6, 7 ]


for title, layers in GENLAYERS:
    layers, txtlayers = build_layers(layers)
    for nproc in GENPROCS:
        for ngll in GENGLL:
            gendict = subst.copy()
            
