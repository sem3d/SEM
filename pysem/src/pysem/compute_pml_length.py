# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to compute PML properties for SEM3D

    Ex.1 : Compute amplitude Ax knowning PML length (300. m):

        python3 compute_pml_length.py @@PML_length 300.

    Ex.2 : Compute PML_length knowing the amplitude Ax (10.)

        python3 compute_pml_length.py @@Ax 10.

Note on SEM3D's material.input "Apow" column
----------------------------------------------
SEM3D's own damping-profile formula (SEM3D/SRC/pml.F90, define_alpha_PML) is
    alpha(x) = Apow * Vp * (1/|L|) * ((x-pos)/L)**npow
i.e. SEM3D already multiplies by the local Vp field and divides by the PML
width L at runtime. The classical CFS-PML target
    sigma_max(L) = -(npow+1)*Vp*ln(R)/(2L)
therefore gives, after solving Apow*Vp/L == sigma_max(L) for Apow, a
dimensionless coefficient that does NOT depend on Vp or L:
    Apow = -(npow+1)*ln(R)/2
This is confirmed by every material.input shipped with SEM3D (SEM3D/TESTS/PMLS/*,
DOC/tutoriel/cas_1): Apow~=10, npow=2 regardless of the material's Vp or the
PML's physical width. Use compute_apow() below for this column instead of
plugging the dimensional sigma_max directly into Apow.
"""
# Required modules
import argparse
import numpy as np

# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, MSSMat UMR CNRS 8579, CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "Cecill-C"
__version__ = "1.0"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"



def compute_apow(npow, RCdb):
    """Dimensionless SEM3D PML amplitude coefficient (material.input "Apow" column).

    Apow = -(npow+1)*ln(R)/2, with R = 10**(RCdb/20) the target reflection
    coefficient. Independent of Vp and PML length: see the module docstring
    for why (SEM3D/SRC/pml.F90 already applies Vp/L at runtime).
    """
    RC = 10.**(RCdb / 20.)
    return -0.5 * (npow + 1) * np.log(RC)


def compute_qmu(vs, coeff=0.1):
    """Shear-wave quality factor from the Qs = coeff*Vs[m/s] rule of thumb."""
    return coeff * vs


def compute_qkappa(vp, vs, qmu, qp_qs_ratio=2.0, min_denom_frac=1e-3):
    """Bulk-modulus quality factor (SEM3D's material.input "Qpression"/Qkappa column).

    Derived like Qlambda but for Kappa = M - (4/3)*mu instead of
    lambda = M - 2*mu (M = rho*Vp**2 the P-wave modulus, mu = rho*Vs**2):
        Qkappa = (Vp**2 - 4/3*Vs**2) / (Vp**2/Qp - 4/3*Vs**2/Qmu)
    Qp is taken as qp_qs_ratio*qmu.

    This weak-attenuation inversion is only physically valid (Qkappa > 0) when
    Vp/Vs > sqrt(4/3*qp_qs_ratio); below that threshold the denominator goes
    negative while the numerator (Kappa, always > 0 for a real solid) stays
    positive, producing a nonsensical negative Qkappa. That threshold is
    1.63 for the default qp_qs_ratio=2.0 -- squarely inside the Vp/Vs~1.15-1.8
    range of ordinary dry consolidated rock, so this is a common case, not a
    rare edge one. When the denominator is negative or too close to zero
    (within min_denom_frac of the numerator scale), fall back to Qkappa=Qp
    (the trivial, always-positive solution when all moduli share one Q) and
    warn instead of returning a negative or blown-up value.
    """
    qp = qp_qs_ratio * qmu
    m = vp**2
    kappa = vp**2 - (4. / 3.) * vs**2
    denom = m / qp - (4. / 3.) * vs**2 / qmu
    if denom <= 0 or denom < min_denom_frac * (m / qp):
        print("Warning: compute_qkappa: Vp/Vs={:.3f} is at or below the "
              "sqrt(4/3*Qp/Qs)={:.3f} stability threshold for Qp/Qs={:.2f}; "
              "the Kappa/Mu split cannot support this Qp,Qmu pair without a "
              "negative Qkappa. Falling back to Qkappa=Qp={:.3f}.".format(
                  vp / vs, np.sqrt(4. / 3. * qp_qs_ratio), qp_qs_ratio, qp))
        return qp
    return kappa / denom


class pml(object):
    def __init__(self,**kwargs):
        self.__call__(**kwargs)
        
    def __call__(self,**kwargs):
        self.__dict__.update(**kwargs)
        self.setup()
        self.check()
    
    def setup(self):
        assert len(self.cp)==len(self.cs)
        self.cp,self.cs = np.array(self.cp),np.array(self.cs)
        self.PML_lengths = np.array(self.PML_lengths)
        self.Lor = self.cs/self.fl[0]                                                       
        self.ll=(self.cs/self.fl[1],self.cp/self.fl[0])
        self.kl=(2.*np.pi/self.ll[1],2.*np.pi/self.ll[0])
        self.RC = 10.**(self.RCdb/20.)
        print("cp: {} m/s".format(self.cp))
        print("cs: {} m/s".format(self.cs))
        print("Frequency limits: {} Hz".format(self.fl))
        print("wave-length limits: {} m".format(self.ll))
        print("wave-number limits: {} 1/m".format(self.kl))
        print("Reflection coefficient: {}".format(self.RC))
    def check(self):
        self.flag = []
        if self.Ax is None:
            self.flag.append('amp')
        if None in self.PML_lengths:
            self.flag.append('len')
        if 'amp' in self.flag:
            if None not in self.PML_lengths:
                self.get_amplitude()
            else:
                raise ValueError('PML lengths not defined!')
        if 'len' in self.flag:
            if self.Ax is not None:
                self.get_length()
            else:
                raise ValueError('PML amplitude not defined!')
            
    def get_length(self):
        if self.PML_type == 'PML':
            self.PML_lengths = (-0.5*((self.p+1)*np.log(self.RC)/self.kl[0]/self.Ax))**(1./(self.p+1))
        elif self.PML_type == 'CPML':
            self.PML_lengths = -0.5/self.Ax*(self.p+1)*self.cp*np.log(self.RC)
        print("PML lengths ({}): {}".format(self.PML_type,self.PML_lengths))
    
    def get_amplitude(self):
        if self.PML_type == 'PML':
            self.Ax = -0.5*(self.p+1)*np.log(self.RC)/(self.kl[0]*self.PML_lengths**(self.p+1))                      
        elif self.PML_type== 'CPML':
            self.Ax = -0.5/self.PML_lengths*(self.p+1)*self.cp*np.log(self.RC)
        print("Ax ({}) = {}".format(self.PML_type,self.Ax))

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@cp',
                        type=float,
                        nargs='+',
                        default=[ 700.,1385.,1732.,3500.],
                        help="P-wave speed in each layer")
    parser.add_argument('@@cs',
                        type=float,
                        nargs='+',
                        default=[ 300., 800.,1000.,2000.],
                        help="S-wave speed in each layer")
    parser.add_argument('@@fl',
                        type=float,
                        nargs='+',
                        default=[0.01,30.],
                        help="Frequency limits")
    parser.add_argument('@@RCdb',
                        type=float,
                        default=-80.,
                        help="Reflection Coefficient in db")
    parser.add_argument('@p',
                        type=int,
                        default=2,
                        help="Polynomial order")
    parser.add_argument('@@Ax',
                        type=float,
                        default=None,
                        help="Ax or d0 coefficient")
    parser.add_argument('@@PML_lengths',
                        type=float,
                        nargs='+',
                        default=None,
                        help="PML length")
    parser.add_argument('@@PML_type',
                        type=str,
                        default="PML",
                        help="PML type [PML|CPML]")
    opt = parser.parse_args().__dict__
    
    p = pml(**opt)    