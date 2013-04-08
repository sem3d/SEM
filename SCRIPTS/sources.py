
from numpy import *
from pylab import *
from matplotlib import mlab
from matplotlib.colors import Normalize

NFFT1=1024
NOVER1=1023
NFFT2=2048
NOVER2=1024

SAVE=True


def Ricker(time, tau, f0):
    sigma = (pi*f0*(time-tau))**2
    return (1-2*sigma)*exp(-sigma)

def Gabor(time, ts, fp, gamma, tau):
    sigma = 2*pi*fp*(time-ts)
    v1 = cos(sigma+pi/2)
    sigma = (sigma/gamma)**2
    v2 = exp(-sigma)
    return v1*v2*tau

def Gaussian(time, ts, tau):
    res = -2*(time-ts)*exp( -(time-ts)**2/tau**2 )
    return res

def Spice_Bench(time, f):
    return (1-(1+time*f)*exp(-time*f))

def plot_specgram(ax, Pxx, f, t, vmin, vmax, tmax, fmax):
    lpxx = 10*log10(Pxx/Pxx.max())
    norm = Normalize(vmin, vmax, clip=True)
    hf_t = (t[1]-t[0])/2.
    hf_f = (f[1]-f[0])/2.
    f = f-hf_f
    t = t-hf_t
    f = concatenate( (f, [f[-1]+hf_f]) )
    t = concatenate( (t, [t[-1]+hf_t]) )
    if fmax:
        print lpxx.shape
        print f.shape
        print t.shape
        fi = f.searchsorted(fmax)
        f = f[:fi+1]
        lpxx = lpxx[:fi,:]
    if tmax:
        ti = t.searchsorted(tmax)
        t = t[:ti+1]
        lpxx = lpxx[:,:ti]
    print fi, ti
    print t
    ax.pcolormesh(t, f, lpxx, norm=norm)

def plot_src(t, s, tmax, fmax, figtitle):
    f=figure(figtitle)
    subplot(3,1,1)
    plot(t,s)
    if tmax: xlim(0, tmax)
    subplot(3,1,2)
    fs = 1./diff(t[:10])[0]
    Pxx, freqs, t = mlab.specgram(s,Fs=fs,NFFT=NFFT1,noverlap=NOVER1)
    plot_specgram(gca(), Pxx, freqs, t, -100., 0. , tmax, fmax)
    colorbar()
    #yscale('log')
    if fmax: ylim(0,fmax)
    #specgram(s,Fs=fs,NFFT=NFFT1,noverlap=NOVER1)
    #colorbar()
    if tmax: xlim(0, tmax)
    subplot(3,1,3)
    psd(s,Fs=fs,NFFT=NFFT2,noverlap=NOVER2)
    if fmax: xlim(0, fmax)
    #gca().set_xscale("log")
    return f

def stat_sig(t,s,name):
    ds = diff(s)
    dt = diff(t)
    print name, " DS[0]=", ds[0], " DT[0]=", dt[0]
    print name, " min/max ds = ", ds.min(), "/", ds.max()

def main():
    t = linspace(0,20.,2**15)
    
    TEST = "TEST0001 : Ricker(0.4,3.)"
    s = Ricker(t, 0.4, 3.)
    plot_src(t, s, 1., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_ricker_1.png")

    TEST = "TEST0002 : Gaussian(0.1, 0.03)"
    s = Gaussian(t, 0.1, 0.03)
    plot_src(t, s, 0.5, 60., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_gaussian.png")
    
    TEST = "TEST0003 : Spice_bench(5.)"
    s = Spice_Bench(t, 2.5)
    plot_src(t, s, 5., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_spice.png")
    figure()
    plot(t[:-1], diff(s))

    s = Ricker(t, 0.2, 4.)
    TEST = "TEST0004 : Ricker(0.2,5.)"
    plot_src(t, s, 1., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_ricker_2.png")
    
    s = Gabor(t, ts=1., fp=3, gamma=2., tau=1.)
    TEST = "Gabor(ts=1,fp=3,g=2,tau=1)"
    plot_src(t, s, 2., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_gabor_1.png")
    
    s = Gabor(t, ts=2., fp=3, gamma=10., tau=1.)
    TEST = "Gabor(ts=4,fp=3,g=10.,tau=1)"
    plot_src(t, s, 4., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_gabor_2.png")
    
    show()

if __name__=="__main__":
    main()
