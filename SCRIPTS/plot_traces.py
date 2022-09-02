#
#
from pylab import *
import h5py
import sys
import argparse


VARIABLES = ["Time", "Displ",
             "Veloc", "Accel", "Pressure", "DUDX" ]

def read_files(files):
    res = []
    for fdesc in files:
        label, capt, fname = fdesc.split(":")
        f = h5py.File(fname)
        data = f[capt][...]
        res.append( (label, capt, data) )
    return res

def plot_dataset(fig, name, capt, data):
    t = data[:,0]
    u = data[:,1:4]
    v = data[:,4:7]
    a = data[:,7:10]
    p = data[:,10]
    dudx = data[:,11]
    ax = fig.add_subplot(3,1,1)
    ax.set_title( "Depla: %s - %s" % (name, capt) )
    ax.plot(t, u[:,0], label="Ux (%s)" % (name,))
    ax.plot(t, u[:,1], label="Uy (%s)" % (name,))
    ax.plot(t, u[:,2], label="Uz (%s)" % (name,))
    ax.legend()
    ax = fig.add_subplot(3,1,2)
    ax.set_title( "Veloc: %s - %s" % (name, capt) )
    ax.plot(t, v[:,0], label="Vx (%s)" % (name,))
    ax.plot(t, v[:,1], label="Vy (%s)" % (name,))
    ax.plot(t, v[:,2], label="Vz (%s)" % (name,))
    ax.legend()
    ax = fig.add_subplot(3,1,3)
    ax.set_title( "Press: %s - %s" % (name, capt) )
    ax.plot(t, dudx[:], label="P (%s)" % (name,))
    ax.legend()


def plot_by_capt(file_descs):
    figs = {}
    for name, capt, data in file_descs:
        fig = figs.get(capt)
        if fig is None:
            fig = figure()
            figs[capt] = fig
        plot_dataset(fig, name, capt, data)

def compare_to_ref(file_descs):
    refs = {}
    for name, capt, data in file_descs:
        if name=="ref":
            refs[capt] = data
    figs = {}
    for name, capt, data in file_descs:
        if name!="ref":
            continue
        fig = figs.get(capt)
        if fig is None:
            fig = figure()
            figs[capt] = fig
        datref = array(data, copy=True)
        datref[:,1:] -= refs[capt][:,1:]
        plot_dataset(fig, name+"-ref", capt, datref)

if __name__=="__main__":
    file_descs = read_files(files)
    plot_by_capt(file_descs)
    compare_to_ref(file_descs)
    savefig("fig.pdf")
    show()
