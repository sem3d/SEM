
from numpy import *
from pylab import *


NFFT1=1024
NOVER1=512
NFFT2=2048
NOVER2=1024

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

def plot_src(t, s, figtitle):
    f=figure(figtitle)
    subplot(3,1,1)
    plot(t,s)
    subplot(3,1,2)
    fs = 1./diff(t[:10])[0]
    specgram(s,Fs=fs,scale_by_freq=True,NFFT=NFFT1,noverlap=NOVER1)
    subplot(3,1,3)
    psd(s,Fs=fs,NFFT=NFFT2,noverlap=NOVER2)
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
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)

    TEST = "TEST0002 : Gaussian(0.1, 0.03)"
    s = Gaussian(t, 0.1, 0.03)
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)
    
    s = Spice_Bench(t, 5.)
    plot_src(t,s, "TEST0003 : Spice_bench(5.)")
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)
    figure()
    plot(t[:-1], diff(s))

    s = Ricker(t, 0.2, 4.)
    TEST = "TEST0004 : Ricker(0.2,5.)"
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)
    
    s = Gabor(t, ts=4., fp=3, gamma=2., tau=1.)
    TEST = "Gabor(4,3,2,1)"
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)
    
    s = Gabor(t, ts=4., fp=3, gamma=10., tau=1.)
    TEST = "Gabor(4,3,10.,1)"
    plot_src(t, s, TEST)
    stat_sig(t, s, TEST)
    
    show()

if __name__=="__main__":
    main()
