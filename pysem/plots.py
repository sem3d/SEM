# -*- coding: utf-8 -*-
# Fonctions de traces pour les capteurs SEM
import os.path as osp
import numpy as np
from pylab import *
from scipy.signal import kaiserord, firwin, lfilter, convolve, medfilt, decimate

class Capteur:
    def __init__(self, name, data, comp):
        self.name = name
        self.data = data
        self.comp = comp
        self.pos = np.array([0.,0.,0.])

    def dist(self):
        return (self.pos**2).sum()

    def set_coord(self, x,y,z):
        self.pos = np.array([x,y,z])

    def __repr__(self):
        return "<{%s %s}>" % (self.name, self.comp)

    def show_stats(self):
        _m = self.data.min(axis=0)
        _M = self.data.max(axis=0)
        _a = self.data.mean(axis=0)
        for comp in 1,2,3:
            print self.name, " XYZ"[comp], ":", _m[comp], " < ", _a[comp], " < ", _M[comp]

    def dt(self):
        # on ne prend pas [1]-[0] car certains filtres decalent le T0
        N=self.data.shape[0]/2
        return self.data[N+1,0] - self.data[N,0]

    def filt_low_pass(self, fmax, fband=2.0):
        sample_rate = 1./self.dt()

        # The Nyquist rate of the signal.
        nyq_rate = sample_rate / 2.0

        # The desired width of the transition from pass to stop,
        # relative to the Nyquist rate.  We'll design the filter
        # with a 5 Hz transition width.
        width = fband/nyq_rate

        # The desired attenuation in the stop band, in dB.
        ripple_db = 60.0

        # Compute the order and Kaiser parameter for the FIR filter.
        N, beta = kaiserord(ripple_db, width)

        # The cutoff frequency of the filter.
        cutoff_hz = fmax

        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        taps = firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))

        # Use lfilter to filter x with the FIR filter.
        data2 = np.array(self.data, copy=True)
        for comp in 1,2,3:
            data2[:,comp] = lfilter(taps, 1.0, self.data[:,comp])

        return self.filt_capteur(data2, "LP(%f Hz)" % fmax)

    def filt_capteur(self, data, pfx):
        cpt = Capteur(self.name+" "+pfx, data, self.comp)
        cpt.pos = self.pos
        return cpt

    def filt_ma(self, n, decimate=True):
        nn = 2*n+1
        flt = np.ones((nn,))/float(nn)
        data = np.array(self.data, copy=True)
        for comp in 1,2,3:
            data[:,comp] = convolve(self.data[:,comp], flt, mode="same")
        if decimate:
            data = np.array(data[n::nn,:], copy=True)
        return self.filt_capteur(data, "MA(%d)" % nn)

    def filt_median(self, n, decimate=True):
        data = np.array(self.data, copy=True)
        nn = 2*n+1
        for comp in 1,2,3:
            data[:,comp] = medfilt(self.data[:,comp], nn)
        if decimate:
            data = np.array(data[n::nn,:], copy=True)
        return self.filt_capteur(data, "MED(%d)" % nn)

    def filt_decimate(self, q):
        data = decimate(self.data, q,axis=0)
        return self.filt_capteur(data, "DEC(%d)" % q)


def parse_capteur_def(f, dir):
    """
    TITRE1
    TITRE2
    NOM_CAPTEUR U0
    FREQ  1
    RAYON 0
    OPERATION  MOY
    GRANDEUR   DEPLA
    COORDX 0.
    COORDY 0.
    COORDZ 0.
    TYPE_CALCUL INTERP
    """
    cpt = None
    name = None
    x = None
    y = None
    z = None
    comp = None
    while True:
        ln = f.readline().strip()
        if ln=="":
            break
        if ln.startswith("NOM_CAPTEUR"):
            name = (ln.split()[1]).strip()
        if ln.startswith("GRANDEUR"):
            comp = (ln.split()[1]).strip().lower()
        if name and comp:
            cpt = read_capteur(name, dir, comp)
        if ln.startswith("COORDX"):
            x = float(ln.split()[1])
        if ln.startswith("COORDY"):
            y = float(ln.split()[1])
        if ln.startswith("COORDZ"):
            z = float(ln.split()[1])
    if x is None or y is None or z is None:
        return None
    cpt.set_coord(x,y,z)
    return cpt


def read_all_capteurs(name="capteurs.dat", dir=None):
    if dir is None:
        dir = "traces"
    f = file(name)
    capteurs = []
    while True:
        cpt = parse_capteur_def(f, dir)
        ln = f.readline()
        if not ln:
            break
        if cpt is None:
            continue
        capteurs.append(cpt)
    return capteurs

def read_depla(name,dir=None):
    """Lecture d'un fichier capteur deplacement
    Le fichier est lu dans traces/ par defaut ou dans dir si présent
    """
    if dir is None:
        dir = "traces"
    return read_capteur(name, dir, "depla")


def read_capteur(name, dir, comp):
    fname = osp.join(dir, name+"_"+comp)
    data = np.fromfile(fname, sep=" ")
    data.shape = -1,4
    return Capteur(name,data,comp)


def plot_capteurs(capts):
    if isinstance(capts,Capteur):
        capts = [capts]
    N = len(capts)
    for comp in 1,2,3:
        for i,c in enumerate(capts):
            subplot(N,3,3*i+comp)
            plot(c.data[:,0], c.data[:,comp], label=c.name+" XYZ"[comp])
            legend()
    #show()

def sel_comp(capts, comp):
    return [c for c in capts if c.comp == comp ]


def stats_capt(cpts):
    for cpt in cpts:
        cpt.show_stats()


def plot_spect(capts, flim=None, detrend=detrend_none, npad=None, NFFT=128):
    if isinstance(capts,Capteur):
        capts = [capts]

    N = len(capts)
    for comp in 1,2,3:
        for i,c in enumerate(capts):
            dt = c.dt()
            Fs = 1./dt
            ax = subplot(2*N,3,6*i+comp)
            plot(c.data[:,0], c.data[:,comp], label=c.name+" XYZ"[comp])
            legend()
            subplot(2*N,3,6*i+comp+3, sharex=ax)
            Pxx, freqs, bins, im = specgram(c.data[:,comp], NFFT=NFFT, Fs=Fs, noverlap=NFFT-1,
                                            detrend=detrend, cmap=cm.gist_heat, pad_to=npad)
            if flim is not None:
                ylim(*flim)
    #show()



def src_spice(t, T):
    return (1-(1+t/T)*np.exp(-t/T))

def src_ricker(time, tau, f0):
    sigma = np.pi * f0 * (time-tau)
    sigma = sigma**2
    return (1 - 2*sigma) * np.exp(-sigma)


