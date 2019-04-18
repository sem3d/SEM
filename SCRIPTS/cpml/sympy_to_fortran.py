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
\\usepackage[margin=1cm]{geometry}
\\usepackage{breqn}
\\begin{document}
"""

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

from sympy import *
#import sympy.galgebra.latex_ex as tex
#tex.Format()
init_printing(use_unicode=True)

a0,a1,a2=symbols("a0 a1 a2", real=True, positive=True)
b0,b1,b2=symbols("beta_0 beta_1 beta_2", real=True, positive=True)
d0,d1,d2=symbols("d0 d1 d2", real=True, positive=True)

#a0,a1,a2=symbols("a0 a1 a2", real=True, positive=True)
#b0,b1,b2=symbols("b0 b1 b2", real=True, positive=True)
#d0,d1,d2=symbols("d0 d1 d2", real=True, positive=True)

w = symbols("omega")
t = symbols("t", real=True, positive=True)

# Calcul des coefficients de la decomposition en elements simples

def s(a,b,w):
    return (b+w)/(a+w)


def get_terms(e,w):
    cst = []
    p = Wild("p")
    for t in e.args:
        for n in (1,2,3):
            u = t.find( 1/(p+w)**n )
            if len(u)==1:
                u = list(u)[0]
                pprint(u)
                print ":->:"
                e = collect(simplify(t/u),(a0,a1,a2))
                pprint(e)
                break
        else:
            cst.append(t)
    print "Const:"
    pprint( Add(*cst))


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

def calc_ijk(lbl, aijk,dijk):
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

    return Lijk_abc, Lijk_aac, Lijk_aba, Lijk_abb, Lijk_aaa
    
abc,aac,aba,abb,aaa = calc_ijk("012", (a0,a1,a2), (d0,d1,d2) )

print "Lijk_abc"
print fcode(abc,source_format="free")
print "Lijk_aac"
print fcode(aac,source_format="free")
print "Lijk_aba"
print fcode(aba,source_format="free")
print "Lijk_abb"
print fcode(abb,source_format="free")
print "Lijk_aaa"
print fcode(aaa,source_format="free")

