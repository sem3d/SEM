# -*- coding: utf-8 -*-
from numpy import *
from pylab import *
from matplotlib import mlab
from matplotlib.colors import Normalize

NFFT1=1024
NOVER1=1023
NFFT2=4096
NOVER2=3072

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

def Spice_Bench(time, f, k=1):
    t = (time*f)**k
    return (1-(1+t)*exp(-t))

def Spice_tanh(time, k, t0):
    return 0.5*(tanh(k*(time-t0)) + 1.)

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
    #colorbar()
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

COUNT=0

def test_ricker(t, tau, f0):
    global COUNT
    TEST = "TEST%02d : Ricker(tau=%f,f0=%f)" % (COUNT, tau, f0)
    s = Ricker(t, tau, f0)
    plot_src(t, s, 1., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_ricker_%02d.png" % COUNT)
    COUNT += 1

def test_gaussian(t, ts, tau):
    global COUNT
    TEST = "TEST%02d : Gaussian(ts=%f, tau=%f)" % (COUNT, ts,tau)
    s = Gaussian(t, ts, tau)
    plot_src(t, s, 0.5, 60., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_gaussiani_%02d.png" % COUNT)
    COUNT += 1

def test_spice(t, f0, gamma):
    global COUNT
    TEST = "TEST%02d : Spice_bench(f0=%f,gamma=%f)" % (COUNT, f0,gamma)
    s = Spice_Bench(t, f0, gamma)
    plot_src(t, s, 5., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_spice%02d.png" %COUNT)
    COUNT += 1

def test_tanh(t, k, t0):
    global COUNT
    TEST = "TEST%02d : Spice_tanh(k=%f,t0=%f)" % (COUNT, k, t0)
    s = Spice_tanh(t, k, t0)
    plot_src(t, s, t0+5., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_tanh%02d.png" % COUNT)
    COUNT += 1

def test_gabor(t, ts, fp, gamma, tau):
    global COUNT
    s = Gabor(t, ts=ts, fp=fp, gamma=gamma, tau=tau)
    TEST = "TEST%02d : Gabor(ts=%f,fp=%f,g=%f,tau=%f)" % (COUNT,ts,fp,gamma,tau)
    plot_src(t, s, 2., 30., TEST)
    stat_sig(t, s, TEST)
    if SAVE: savefig("src_gabor_%02d.png" % COUNT)
    COUNT += 1


def main():
    t = linspace(0,20.,2**15)
    #test_ricker(t, 0.4, 3.)
    #test_gaussian(t, 0.1, 0.03)
    #test_spice(t, 2.5, 1.)
    #test_spice(t, 0.75, 1.)
    #test_ricker(t, 0.2, 4.)
    #test_gabor(t, ts=1., fp=3., gamma=2., tau=1.)
    #test_gabor(t, ts=2., fp=3., gamma=10., tau=1.)
    test_spice(t, 0.75, 2.)
    test_tanh(t, 4, 5.)

    show()

if __name__=="__main__":
    main()
