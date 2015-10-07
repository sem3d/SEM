# -*- mode: python; coding: utf-8 -*-
# Outil pour la génération de cas de non-régression

import sys
import os
import os.path as osp

template = sys.argv[1]
outputdir = sys.argv[2]

MATERIAL = [
    ("S", 6300., 2500., 2800., 5, 1000., 1000.), # S/F Vp Vs Rho ngll Qp Qs
    ("F", 1300.,    0., 1000., 5, 1000., 1000.),
]

MATDEF="%(t)s %(vp)s %(vs)s %(rho)s  %(ngll)s  %(ngll)s  %(ngll)s 0.1 %(qp)s %(qs)s"
TXTLAYER="Epaisseur %(laythick)s m %(ndiv)s mailles"
LAYER="%(laythick)s %(ndiv)s # Layer %(laynum)s"

MAT_DAT_TEMPLATE = open(osp.join(template,"mat.dat")).read()
MATER_IN_TEMPLATE = open(osp.join(template,"mater.in")).read()
INPUT_SPEC_TEMPLATE = open(osp.join(template,"input.spec")).read()
CAPTEURS_DAT_TEMPLATE = open(osp.join(template,"capteurs.dat")).read()

subst = {
    "name" : "Test_name",
    "save_snap" : "true",
    "snap_interval" : 0.01,
    "prorep" : "false",
    "restart_iter" : 0,
    "sim_time" : 0.5,

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
# Attenuation
    "nsolids" : 0,
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
    ("Solid/noATN/noPML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                    "subst" : {"pmla" : 0, "nsolids":0},   } ),
    ("Fluid/noATN/noPML/g%(nglls)s/case", { "layers": [ (600, 12, MATERIAL[1]) ],
                                    "subst" : {"pmla" : 0, "nsolids":0},}),
    ("SolidFluid/noATN/noPML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                    "subst" : {"pmla" : 0, "nsolids":0},}),

    ("Solid/noATN/PML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                  "subst" : {"pmla" : 1, "nsolids":0},   } ),
    ("Fluid/noATN/PML/g%(nglls)s/case", { "layers": [ (600, 12, MATERIAL[1]) ],
                                  "subst" : {"pmla" : 1, "nsolids":0},}),
    ("SolidFluid/noATN/PML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                  "subst" : {"pmla" : 1, "nsolids":0},}),

    ("Solid/ATN1/noPML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                    "subst" : {"pmla" : 0, "nsolids":1},   } ),
    ("SolidFluid/ATN1/noPML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                    "subst" : {"pmla" : 0, "nsolids":1},}),

    ("Solid/ATN1/PML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                 "subst" : {"pmla" : 1, "nsolids":1},   } ),
    ("SolidFluid/ATN1/PML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                 "subst" : {"pmla" : 1, "nsolids":1},}),

    ("Solid/ATN3/noPML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                    "subst" : {"pmla" : 0, "nsolids":3},   } ),
    ("SolidFluid/ATN3/noPML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                    "subst" : {"pmla" : 0, "nsolids":3},}),

    ("Solid/ATN3/PML/g%(nglls)s/case", { "layers" : [ (600, 12, MATERIAL[0]) ],
                                 "subst" : {"pmla" : 1, "nsolids":3},   } ),
    ("SolidFluid/ATN3/PML/g%(nglls)s/case", { "layers": [ (300,  6, MATERIAL[0]), (300, 6, MATERIAL[1]) ],
                                 "subst" : {"pmla" : 1, "nsolids":3},}),
]


class SemCase(object):
    def __init__(self, casename, template, casespec):
        self.casename = casename
        self.template = template.copy()
        self.parse_casespec(casespec)

    def parse_casespec(self, casespec):
        self.compute_layers(casespec)
        self.compute_material(casespec)

    def compute_layers(self, casespec):
        layers = []
        for k, layer_def in enumerate(casespec):
            dct = { "laythick" : layer_def[0], "ndiv": layer_def[1], "laynum":k }
            layers.append( LAYER % dct )
        self.template["nblayers"] = len(layers)
        self.template["layers"] = "\n".join(layers)

    def compute_material(self, casespec):
        layers = []
        for k, layer_def in enumerate(casespec):
            ltype, vp, vs, rho, ngll, qp, qs = layer_def[2]
            dct = { "t" : ltype,
                    "vp" : vp,
                    "vs":vs,
                    "rho":rho,
                    "ngll" : ngll,
                    "qp" : qp,
                    "qs" : qs,
                    }
            layers.append( MATDEF % dct )
        self.template["nmat"] = len(layers)
        self.template["matdef"] = "\n".join(layers)
        
    def gen_case_dir(self, destdir):
        dirname = osp.join(destdir, self.casename)
        try:
            os.makedirs(dirname)
        except OSError, err:
            if err.errno!=17:
                raise
        self.gen_mater_in(dirname)
        self.gen_mat_dat(dirname)
        self.gen_input_spec(dirname)
        self.gen_capteurs_dat(dirname)


    def gen_mater_in(self, dirname):
        dest = open(osp.join(dirname, "mater.in"), "w")
        dest.write(MATER_IN_TEMPLATE % self.template)
        

    def gen_mat_dat(self, dirname):
        dest = open(osp.join(dirname, "mat.dat"), "w")
        dest.write(MAT_DAT_TEMPLATE % self.template)

    def gen_input_spec(self, dirname):
        dest = open(osp.join(dirname, "input.spec"), "w")
        dest.write(INPUT_SPEC_TEMPLATE % self.template)


    def gen_capteurs_dat(self, dirname):
        dest = open(osp.join(dirname, "capteurs.dat"), "w")
        dest.write(CAPTEURS_DAT_TEMPLATE % self.template)



def main():
    subst["nglls"] = 5
    for testname, testdesc in GENLAYERS:
        name = testname % subst
        case = SemCase(name, subst, testdesc["layers"])
        case.gen_case_dir(outputdir)

if __name__ == "__main__":
    main()
