# Affiche des info de debugs sur les communications SEM
import h5py
import sys
import os.path as osp
import pylab as plt
import numpy as np

nprocs = int(sys.argv[1])

NGLL=5

def show_comm_info(n, fname):
    f = h5py.File(fname,"r")
    ncomms = f.attrs["tot_comm_proc"]
    sizes = []
    sizetot = 0
    if n%16==0:
        print "Proc: %5d" %n
    for nc in range(ncomms):
        grp = f["Proc%04d" % nc]
        nfaces = grp.attrs["n_faces"]
        nedges = grp.attrs["n_edges"]
        nverts = grp.attrs["n_vertices"]
        dest = grp.attrs["proc_dest"]
        nbytes = 3*8*( nfaces*(NGLL-2)*(NGLL-2) + nedges*(NGLL-2) + nverts)
        sizes.append(nbytes)
        sizetot += nbytes
        #print "-> %5d : %4d/%4d/%4d f/e/v : ~ %8d bytes" % (dest, nfaces, nedges, nverts, nbytes)
    f.close()
    return sizes, sizetot, ncomms

sizes = []
procsizes = []
nc = []
for k in range(nprocs):
    fname = osp.join("sem", "mesh4spec.%04d.h5"%k)
    psizes, psz, pnc =show_comm_info(k, fname)
    sizes += psizes
    procsizes.append(psz)
    nc.append(pnc)

sizes = np.array(sizes)
nc = np.array(nc)
procsizes=np.array(procsizes)


print "Messages: min/max", sizes.min(), sizes.max()
print "PerProc: min/max", procsizes.min(), procsizes.max()
print "Comms:", nc.min(), nc.max()



plt.subplot(3,1,1)
plt.hist(sizes, bins=max(10, nprocs/20))
plt.subplot(3,1,2)
plt.hist(procsizes, bins=max(10, nprocs/20))
plt.subplot(3,1,3)
plt.plot(nc)
plt.show()
