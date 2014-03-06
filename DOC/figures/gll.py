
from numpy.linalg import inv

GLL5x=[ -1., -0.65465367070797720, 0., 0.65465367070797720, 1.0   ]

GLL5w=[ 0.10000000000000001, 0.54444444444444440, 0.71111111111111114, 0.54444444444444440, 0.10000000000000001]   

  
from pylab import *
GLL9x=[ -1.0000000000000000,   -0.89975799541146018,   -0.67718627951073773,   -0.36311746382617816, 0.0000000000000000,    0.36311746382617816,    0.67718627951073773,    0.89975799541146018, 1.0000000000000000,  ]

GLL9w=[  2.7777777777777776E-002,  0.16549536156080558,    0.27453871250016165,    0.34642851097304639,    0.37151927437641719,    0.34642851097304639,    0.27453871250016165,    0.16549536156080558, 2.7777777777777776E-002]


X,Y = meshgrid(GLL9x,GLL9x)
figure()
for x in GLL9x:
    axhline(x, -1.,1.)
    axvline(x, -1.,1.)
    plot(X.flat, Y.flat, "ro")

tight_layout()
savefig("elem_gll_2d.png")

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



t=linspace(-1,1,500)

figure()
plag=[]
for i in range(9):
    p = pol(t,i,GLL9x)
    plot(t,p)
    plag.append(p)

tight_layout()
savefig("lagrange_9.png")

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
savefig("shape_fct_9x9.png")


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
savefig("mass_9x9.png")


show()



