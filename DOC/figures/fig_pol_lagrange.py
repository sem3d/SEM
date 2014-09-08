
from numpy.linalg import inv
from matplotlib import pyplot as plt
import matplotlib.font_manager
matplotlib.font_manager.USE_FONTCONFIG=True

plt.xkcd()

GLL5x=[ -1., -0.65465367070797720, 0., 0.65465367070797720, 1.0   ]

GLL5w=[ 0.10000000000000001, 0.54444444444444440, 0.71111111111111114, 0.54444444444444440, 0.10000000000000001]

from pylab import *
GLL9x=[ -1.0000000000000000,   -0.89975799541146018,   -0.67718627951073773,   -0.36311746382617816, 0.0000000000000000,    0.36311746382617816,    0.67718627951073773,    0.89975799541146018, 1.0000000000000000,  ]

GLL9w=[  2.7777777777777776E-002,  0.16549536156080558,    0.27453871250016165,    0.34642851097304639,    0.37151927437641719,    0.34642851097304639,    0.27453871250016165,    0.16549536156080558, 2.7777777777777776E-002]


def figure_elem_gll(figname):
    X,Y = meshgrid(GLL9x,GLL9x)
    figure()
    for x in GLL9x:
        axhline(x, -1.,1.)
        axvline(x, -1.,1.)
        plot(X.flat, Y.flat, "ro")

    tight_layout()
    savefig(figname)

def pol(x, i, cx):
    p = 1.
    q = 1.
    n = len(cx)
    for j in range(n):
        if i==j:
            continue
        p = p*(x-cx[j])
        q = q*(cx[i]-cx[j])
    return p/q

def all_poly(t, glls):
    plag=[]
    for i in range(len(glls)):
        p = pol(t,i,glls)
        plag.append(p)
    return plag

def figure_poly_lagrange(t, plag, figname):
    N = len(plag)
    figure()
    for i in range(N):
        plot(t,plag[i])

    tight_layout()
    savefig(figname)

def figure_2d_shape_func(plag, figname):
    figure()
    idx=[1,2,4]

    k=1
    for i in idx:
        for j in idx:
            subplot(3,3,k)
            k=k+1
            qi = array(plag[i],copy=True)
            qj = array(plag[j],copy=True)
            qi.shape = 1,-1
            qj.shape = -1,1
            P = qi*qj
            imshow(P,extent=[-1,1,-1,1],vmin=-1,vmax=1.)
            title("$\Phi_%d(x)*\Phi_%d(y)$"%(i,j))

    tight_layout()
    savefig(figname)


def figure_mass_mat(plag, figname):
    figure()

    idx=range(9)
    M = np.zeros( (9,9), float)
    for i in idx:
        for j in idx:
            qi = array(plag[i],copy=True)
            qj = array(plag[j],copy=True)
            M[i,j] = (qi*qj).sum()

    mm = abs(M).max()
    pcolor(M,vmin=-mm, vmax=mm)
    colorbar()
    title("$M_{i,j}$")
    figure()
    N = inv(M)
    mm = abs(N).max()
    pcolor(N,vmin=-mm, vmax=mm)
    colorbar()

    tight_layout()
    savefig(figname)




t=linspace(-1,1,500)
plag9 = all_poly(t, GLL9x)
plag5 = all_poly(t, GLL5x)
figure_elem_gll("elem_gll_2d.pdf")
figure_poly_lagrange(t, plag9, "lagrange_9.pdf")
figure_poly_lagrange(t, plag5, "lagrange_5.pdf")
figure_2d_shape_func(plag9, "shape_fct_9x9.pdf")
figure_mass_mat(plag9, "mass_9x9.pdf")
show()
