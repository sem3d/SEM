# -*- coding: utf-8 -*-
import sys
import os
import h5py
import optparse

PATTERN = "mesh4spec.%04d.h5"

def find_mesh_files(pattern):
    files = []
    proc = 0
    while True:
        fname = pattern % proc
        print "Trying:", fname
        if not os.path.exists(fname):
            break
        files.append((proc,fname))
        proc = proc+1
    return files

XMF_PREC="""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
<Grid CollectionType="Spatial" GridType="Collection" Name="Proc">
"""

def write_xmf_file(output, files):
    xmf = open(output, "w")
    xmf.write(XMF_PREC)

    for nproc, fname in files:
        f = h5py.File(fname)
        NE = f["elements"].shape[0]
        NN = f["local_nodes"].shape[0]

        xmf.write('<Grid Name="%s">\n' % fname)
        xmf.write('<Topology Type="Hexahedron" NumberOfElements="%d">\n' % NE)
        xmf.write('<DataItem Format="HDF" DataType="Int" Dimensions="%d 8">' % NE)
        xmf.write('%s:/elements\n</DataItem>\n' % fname)
        xmf.write('</Topology>\n')
        xmf.write('<Geometry Type="XYZ">\n')
        xmf.write('<DataItem Format="HDF" DataType="Int" Dimensions="%d 3">' % NN)
        xmf.write('%s:/local_nodes\n</DataItem>\n' % fname)
        xmf.write('</Geometry>\n')

        xmf.write('<Attribute Name="Mat" Center="Cell" AttributeType="Vector" Dimensions="%d 3">\n' % NE)
        xmf.write('<DataItem Format="HDF" DataType="Int" Dimensions="%d 3">' % NE)
        xmf.write('%s:/material\n</DataItem>\n' % fname)
        xmf.write('</Attribute>\n')

        xmf.write('</Grid>\n')
    xmf.write('</Grid>\n</Domain>\n</Xdmf>\n')


def main():
    parser = optparse.OptionParser()
    parser.add_option("-p", dest="pattern", default=PATTERN, action="store",
                      help="Pattern for mesh2spec files")
    opts, args=parser.parse_args()
    output = args[0]

    files  = find_mesh_files(opts.pattern)
    write_xmf_file(output, files)


if __name__=="__main__":
    main()
