import sys
import h5py

NP = int(sys.argv[1])


for n in range(NP):
    f = h5py.File("sem/mesh4spec.%04d.h5"%n)
    nc = f.attrs["tot_comm_proc"]
    lst = []
    for k in range(nc):
        dst = f["Proc%04d"%k].attrs["proc_dest"]
        lst.append(str(dst))
    print n,":",",".join(lst)
