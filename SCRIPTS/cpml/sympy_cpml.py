# -*- coding: utf-8 -*-

################ Readline stuff for interactive sessions
import os
import atexit

histfile = os.path.join(os.environ["HOME"],".pyhist")

try:
    import readline
    try:
        readline.read_history_file(histfile)
    except IOError:
        pass
    atexit.register(readline.write_history_file, histfile)
except ImportError:
    print "Readline not available"
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")

del os, histfile
#####################

HEADER_TEX = """
\\documentclass[a4paper,10pt]{article}
\\usepackage[margin=1cm,landscape]{geometry}
%\\usepackage{breqn}
\\usepackage{amsmath}
\\begin{document}
"""
def write_as_terms(f,e,w):
    cst = []
    d = {}
    p = Wild("p")
    for t in e.args:
        for n in (1,2,3):
            u = t.find( 1/(p+w)**n )
            if len(u)==1:
                u = list(u)[0]
                d[u] = collect(simplify(t/u),(a0,a1,a2))
                break
        else:
            cst.append(t)
    cterm = Add(*cst)
    f.write(" &=" + latex(cterm))
    # Trick to sort terms in w
    if len(d)>1:
        keys = Add(*d.keys()).args
    else:
        # With one Item Add(x) reduces to x and x.args is not what we want
        keys = d.keys()
    for u in keys:
        f.write("\\\\\n")
        f.write("&+ ")
        f.write(latex(d[u]))
        f.write(" \\cdot ")
        f.write(latex(u))

def write_def(f, lbl, expr):
    f.write("\n\\begin{eqnarray}\n %s & = & " % lbl)
    terms = []
    tail = expr
    while isinstance(tail,Add):
        head, tail = tail.as_two_terms()
        terms.append(head)
    for t in terms:
        f.write(latex(t))
        f.write("\\\\\n & + & ")
    f.write(latex(tail))
    f.write("\n\\end{eqnarray}\n")

def write_def(f, lbl, expr):
    f.write("\n\\begin{dmath*}\n %s = " % lbl)
    f.write(latex(expr))
    f.write("\n\\end{dmath*}\n")

def write_def(f, lbl, expr):
    f.write("\n\\begin{align*}\n %s " % lbl)
    write_as_terms(f, expr, w)
    f.write("\n\\end{align*}\n")

from sympy import *
#import sympy.galgebra.latex_ex as tex
#tex.Format()
init_printing(use_unicode=True)

a0,a1,a2=symbols("alpha_0 alpha_1 alpha_2", real=True, positive=True)
b0,b1,b2=symbols("beta_0 beta_1 beta_2", real=True, positive=True)
d0,d1,d2=symbols("delta_0 delta_1 delta_2", real=True, positive=True)

#a0,a1,a2=symbols("a0 a1 a2", real=True, positive=True)
#b0,b1,b2=symbols("b0 b1 b2", real=True, positive=True)
#d0,d1,d2=symbols("d0 d1 d2", real=True, positive=True)

w = symbols("omega")
t = symbols("t", real=True, positive=True)

# Calcul des coefficients de la decomposition en elements simples

def s(a,b,w):
    return (b+w)/(a+w)

# Calcul de L = F-1(w**2.s0.s1.s2)  (EA)

# 3 racines distinces
eL_abc = w*w*s(a0,a0+d0,w)*s(a1,a1+d1,w)*s(a2,a2+d2,w)
L_abc = apart(eL_abc,w)

# a0==a1
eL_aac = w*w*s(a0,a0+d0,w)*s(a0,a0+d1,w)*s(a2,a2+d2,w)
L_aac = apart(eL_aac,w)

# a0==a2
eL_aba = w*w*s(a0,a0+d0,w)*s(a1,a1+d1,w)*s(a0,a0+d2,w)
L_aba = apart(eL_aba,w)

# a1==a2
eL_abb = w*w*s(a0,a0+d0,w)*s(a1,a1+d1,w)*s(a1,a1+d2,w)
L_abb = apart(eL_abb,w)

# a0==a1==a2
eL_aaa = w*w*s(a0,a0+d0,w)*s(a0,a0+d1,w)*s(a0,a0+d2,w)
L_aaa = apart(eL_aaa,w)


f = file("output.tex","w")
f.write(HEADER_TEX)

f.write("\\section{Calcul des termes en L pour Solide et Fluide}\n")

# 3 racines distinces
eL_ab = w*w*s(a0,a0+d0,w)*s(a1,a1+d1,w)
L_ab = apart(eL_abc,w)

# a0==a1
eL_aa = w*w*s(a0,a0+d0,w)*s(a0,a0+d1,w)
L_aa = apart(eL_aa,w)

f.write("\\subsection{$L = \\mathcal{F}^{-1}(s_{0} s_{1} s_{2})$}\n")
write_def(f, "L_{abc}", L_abc)
write_def(f, "L_{aac}", L_aac)
write_def(f, "L_{aba}", L_aba)
write_def(f, "L_{abb}", L_abb)
write_def(f, "L_{aaa}", L_aaa)

f.write("\\subsection{$L = \\mathcal{F}^{-1}(s_{0} s_{1})$}\n")
write_def(f, "L_{ab}", L_ab)
write_def(f, "L_{aa}", L_aa)

def calc_ijk(f, lbl, aijk,dijk):
    ai,aj,ak = aijk
    di,dj,dk = dijk

    # Calcul de Lijk = F-1(si.sj/sk)  (Lijk)

    # 3 racines distinctes
    eLijk_abc = s(ai,ai+di,w)*s(aj,aj+dj,w)/s(ak,ak+dk,w)
    Lijk_abc = apart(eLijk_abc,w)

    # a0==a1
    eLijk_aac = s(ai,ai+di,w)*s(ai,ai+dj,w)/s(ak,ak+dk,w)
    Lijk_aac = apart(eLijk_aac,w)

    # a0==b2
    eLijk_aba = s(ai,ai+di,w)*s(aj,aj+dj,w)/s(ak,ai,w)
    Lijk_aba = apart(eLijk_aba,w)

    # a1==b2
    eLijk_abb = s(ai,ai+di,w)*s(aj,aj+dj,w)/s(ak,aj,w)
    Lijk_abb = apart(eLijk_abb,w)

    # a0==a1==b2
    eLijk_aaa = s(ai,ai+di,w)*s(ai,ai+dj,w)/s(ak,ai,w)
    Lijk_aaa = apart(eLijk_aaa,w)

    f.write("\\subsection{$L^{%s} = \\mathcal{F}^{-1}(\\frac{s_{%s}s_{%s}}{s_{%s}})$}\n\n" % (lbl, lbl[0],lbl[1],lbl[2]))
    write_def(f, "L^{%s}_{abc}" % lbl, Lijk_abc)
    write_def(f, "L^{%s}_{aac}" % lbl, Lijk_aac)
    write_def(f, "L^{%s}_{aba}" % lbl, Lijk_aba)
    write_def(f, "L^{%s}_{abb}" % lbl, Lijk_abb)
    write_def(f, "L^{%s}_{aaa}" % lbl, Lijk_aaa)


f.write("\\section{Calcul des termes $L_{ijk}$}\n")

calc_ijk(f, "012", (a0,a1,a2), (d0,d1,d2) )
calc_ijk(f, "021", (a0,a2,a1), (d0,d2,d1) )
calc_ijk(f, "120", (a1,a2,a0), (d1,d2,d0) )


def calc_ixj(f, lbl, aij,dij):
    ai,aj = aij
    di,dj = dij

    # Calcul de Lijk = F-1(si.sj/(sk=1))  (Lijk)

    # 3 racines distinctes
    eLijk_abc = s(ai,ai+di,w)*s(aj,aj+dj,w)
    Lijk_abc = apart(eLijk_abc,w)

    # a0==a1
    eLijk_aac = s(ai,ai+di,w)*s(ai,ai+dj,w)
    Lijk_aac = apart(eLijk_aac,w)

    f.write("\\subsection{$L^{%s} = \\mathcal{F}^{-1}(s_{%s}s_{%s})$}\n\n" % (lbl, lbl[0],lbl[1]))
    write_def(f, "L^{%s}_{ab}" % lbl, Lijk_abc)
    write_def(f, "L^{%s}_{aa}" % lbl, Lijk_aac)

def calc_ioj(f, lbl, aij, dij):
    ai,aj = aij
    di,dj = dij

    # Calcul de Lijk = F-1(si.1/sj)  (Lijk)

    # 3 racines distinctes
    eLijk_abc = s(ai,ai+di,w)/s(aj,aj+dj,w)
    Lijk_abc = apart(eLijk_abc,w)

    # a0==a1
    eLijk_aac = s(ai,ai+di,w)/s(aj,ai,w)
    Lijk_aac = apart(eLijk_aac,w)

    f.write("\\subsection{$L^{%s} = \\mathcal{F}^{-1}(\\frac{s_{%s}}{s_{%s}})$}\n\n" % (lbl, lbl[0],lbl[1]))
    write_def(f, "L^{%s}_{ab}" % lbl, Lijk_abc)
    write_def(f, "L^{%s}_{aa}" % lbl, Lijk_aac)

calc_ixj(f, "01x", (a0,a1), (d0,d1) )
calc_ixj(f, "02x", (a0,a2), (d0,d2) )
calc_ixj(f, "12x", (a1,a2), (d1,d2) )

calc_ioj(f, "0x2", (a0,a2), (d0,d2) )
calc_ioj(f, "0x1", (a0,a1), (d0,d1) )
calc_ioj(f, "1x0", (a1,a0), (d1,d0) )




f.write("\\section{Calcul des termes en L pour le couplage}\n")

# On traite F-1(s0s1s2/si) sachant qu'on a au plus s0 et s1
f.write("\\subsection{Terme : $\\mathcal{F}^{-1}(s_{0})$}\n")

eL0 = s(a0,a0+d0,w)
L0 = apart(eL0,w)
write_def(f, "L^{0}", L0)

f.write("\\subsection{Terme  ab : $\\mathcal{F}^{-1}(s_{0}s_{1})$}\n")
eL01_ab = s(a0,a0+d0,w)*s(a1,a1+d1,w)
L01_ab = apart(eL01_ab,w)
write_def(f, "L^{01}_{ab}", L01_ab)

f.write("\\subsection{Terme  aa : $\\mathcal{F}^{-1}(s_{0}s_{1})$}\n")
eL01_aa = s(a0,a0+d0,w)*s(a0,a0+d1,w)
L01_aa = apart(eL01_aa,w)
write_def(f, "L^{01}_{aa}", L01_aa)

# On traite F-1(s0s1s2/si) sachant qu'on a au plus s0 et s1
f.write("\\subsection{Terme : $\\mathcal{F}^{-1}(\omega^{2}s_{0})$}\n")

eL0 = w*w*s(a0,a0+d0,w)
L0 = apart(eL0,w)
write_def(f, "L^{0}", L0)

f.write("\\subsection{Terme  ab : $\\mathcal{F}^{-1}(\omega^{2}s_{0}s_{1})$}\n")
eL01_ab = w*w*s(a0,a0+d0,w)*s(a1,a1+d1,w)
L01_ab = apart(eL01_ab,w)
write_def(f, "L^{01}_{ab}", L01_ab)

f.write("\\subsection{Terme  aa : $\\mathcal{F}^{-1}(\omega^{2}s_{0}s_{1})$}\n")
eL01_aa = w*w*s(a0,a0+d0,w)*s(a0,a0+d1,w)
L01_aa = apart(eL01_aa,w)
write_def(f, "L^{01}_{aa}", L01_aa)


f.write("\\end{document}\n")
#tex.xdvi()
