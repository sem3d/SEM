# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to create h5 material files for SEM3D

    Ex.1 : Compute amplitude Ax knowning PML length (300. m):
        
        python3 generate_h5_materials.py --pfx base_smooth_het
"""
# Required modules
import argparse
import h5py
import numpy as np

# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

dirs_dict = {'x':0,'y':1,'z':2}
trnsp = (2,1,0)

def base_smooth_heterogeneous(d,grd,nu=0.3):
    z = grd[d]
    la = np.zeros_like(grd[d])
    mu = np.zeros_like(grd[d])
    ds = np.full_like(grd[d],2000.).transpose(*trnsp)
    
    la = 20.*(80.0+0.45*np.abs(z)+\
        35.0*np.exp(-(np.abs(z)-22.5)**2/150.0))*1.e6 # Mpa
    mu = 0.5*(1.-2.*nu)*la/nu
    la = la.transpose(*trnsp)
    mu = mu.transpose(*trnsp)
    vp = np.sqrt((la+2.*mu)/ds)
    vs = np.sqrt(mu/ds)
    return {'la':la,'mu':mu,'ds':ds,'vp':vp,'vs':vs}

def linear_gradient(d,grd,nu=0.3):
    z = grd[d]
    la = np.zeros_like(grd[d])
    mu = np.zeros_like(grd[d])
    ds = np.full_like(grd[d],2000.).transpose(*trnsp)
    
    la = (100.0+0.45*np.abs(z))*1.e6 # Mpa
    mu = 0.5*(1.-2.*nu)*la/nu
    la = la.transpose(*trnsp)
    mu = mu.transpose(*trnsp)
    vp = np.sqrt((la+2.*mu)/ds)
    vs = np.sqrt(mu/ds)
    return {'la':la,'mu':mu,'ds':ds,'vp':vp,'vs':vs}

func = {'base_smooth_heterogeneous':base_smooth_heterogeneous,
        'linear_gradient':linear_gradient}

def grid(lims):
    xv = np.linspace(lims['xmin'],lims['xmax'],lims['nx'],dtype=np.float64)
    yv = np.linspace(lims['ymin'],lims['ymax'],lims['ny'],dtype=np.float64)
    zv = np.linspace(lims['zmin'],lims['zmax'],lims['nz'],dtype=np.float64)
    xg,yg,zg = np.meshgrid(xv,yv,zv,indexing='xy')
    return (xg,yg,zg)

def gen_mat(model,dirs,prop,grd,nu=0.3):
    d = dirs_dict[dirs]
    mats = func[model.lower()](d,grd,nu)
    return dict(tuple([(v,mats[v]) for v in prop]))

def write_h5(pfx,prop,mat,lims,xdmf=True):
    for v in prop:
        with h5py.File("{}_{}.h5".format(pfx,v),"w") as fid:
            samples = fid.create_dataset("samples",mat[v].shape,data=mat[v])
            # chunks=tuple([l//10 for l in mat[v].shape]))
            xMinGlob = np.array([lims['xmin'],lims['ymin'],lims['zmin']])
            xMaxGlob = np.array([lims['xmax'],lims['ymax'],lims['zmax']])
            fid.attrs['xMinGlob'] = xMinGlob
            fid.attrs['xMaxGlob'] = xMaxGlob
            fid.close()

def write_xdmf(pfx,prop,mat,lims):
    
    xMinGlob = np.array([lims['xmin'],lims['ymin'],lims['zmin']])
    xMaxGlob = np.array([lims['xmax'],lims['ymax'],lims['zmax']])
    dxV = (xMaxGlob-xMinGlob)/np.array([lims['nx']-1,lims['ny']-1,lims['nz']-1])
                               
    for v in prop:
        szs = mat[v].shape
        with open("{}_{}.xmf".format(pfx,v),"w") as fid:
            fnm = "{}_{}".format(pfx,v)
            fid.write('''<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n'''+
                      '''<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n'''+
                      '''<Domain>\n''')
            fid.write('   <DataItem Name="{}" Format="HDF" DataType="Float" Precision="8" Dimensions="{} {} {}">\n'.format(v,*szs))
            fid.write('        %s.h5:/samples\n'%fnm)
            fid.write('   </DataItem>\n')
            fid.write('  <Grid GridType="Collection" CollectionType="Spatial">\n')
            fid.write('   <Grid Name="Group1">\n')
            fid.write('     <Topology TopologyType="3DCoRectMesh" Dimensions="%u %u %u"/>\n'%szs)
            fid.write('     <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
            fid.write('   <DataItem Name="origin" Format="XML" DataType="Float" Precision="8" Dimensions="3">\n')
            fid.write('         %30.10f\n'%lims['zmin'])
            fid.write('         %30.10f\n'%lims['ymin'])
            fid.write('         %30.10f\n'%lims['xmin'])
            fid.write('   </DataItem>\n')
            fid.write('   <DataItem Name="step" Format="XML" DataType="Float" Precision="8" Dimensions="3">\n')
            fid.write('         %30.10f\n'%dxV[2])
            fid.write('         %30.10f\n'%dxV[1])
            fid.write('         %30.10f\n'%dxV[0])
            fid.write('   </DataItem>\n')
            fid.write('     </Geometry>\n')
            fid.write('     <Attribute Name="%s" Center="Node" AttributeType="Scalar">\n'%v)
            fid.write('       <DataItem Reference="XML">\n')
            fid.write('         /Xdmf/Domain/DataItem[@Name="%s"]\n'%v)
            fid.write('       </DataItem>\n')
            fid.write('     </Attribute>\n')
            fid.write('   </Grid>\n')
            fid.write('  </Grid>\n')
            fid.write(' </Domain>\n')
            fid.write('</Xdmf>\n')
            fid.close()

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--prop',type=str,nargs='+',default= ['la','mu','ds','vp','vs'],help="list of properties to be generated")
    parser.add_argument('--tag',type=str,default="linear_gradient",help="tag for material model")
    parser.add_argument('--dir',type=str,default="z",help="Main gradient direction [x|y|z]")
    parser.add_argument('--xlim',type=float,nargs='*',default=[-45.0,45.0],help="Limits of the box [xmin xmax]")
    parser.add_argument('--ylim',type=float,nargs='*',default=[-45.0,45.0],help="Limits of the box [ymin ymax]")
    parser.add_argument('--zlim',type=float,nargs='*',default=[-65.0, 5.0],help="Limits of the box [zmin zmax]")
    parser.add_argument('--step',type=float,nargs='*',default=[10,10,501],help="Numbers of points per direction [nx ny nz]")
    parser.add_argument('--pfx',type=str,default="linear_gradient",help="File prefix")
    parser.add_argument('--nu',type=float,default=0.3,help="Poisson's ratio")
    opt = parser.parse_args().__dict__
    
    assert len(opt['xlim'])==2
    assert len(opt['ylim'])==2
    assert len(opt['zlim'])==2
    
    opt['xlim'].sort() 
    opt['ylim'].sort()
    opt['zlim'].sort()

    # mechanical properties to be generated
    model = (opt['tag'],opt['dir'])
    print("Model tag: {0} - Model grad: {1}".format(*model))
    
    # define grid limits (larger than SEM3D domain)
    lims = (('xmin',opt['xlim'][0]),('xmax',opt['xlim'][1]),
            ('ymin',opt['ylim'][0]),('ymax',opt['ylim'][1]),
            ('zmin',opt['zlim'][0]),('zmax',opt['zlim'][1]),
            ('nx',opt['step'][0]),('ny',opt['step'][1]),('nz',opt['step'][2]))
    lims = dict(lims)
    print(f'Domain limits/discretization : \n X: {lims["xmin"]} m : {lims["xmax"]} m; nx={lims["nx"]} \n Y: {lims["ymin"]} m : {lims["ymax"]} m; ny={lims["ny"]} \n Z: {lims["zmin"]} m : {lims["zmax"]} m; nz={lims["nz"]}')
    
    # generate grid points
    grd = grid(lims)
    
    # generate model
    mat = gen_mat(*model,opt['prop'],grd,opt['nu'])
    
    # write hdf5 file
    write_h5(opt['pfx'],opt['prop'],mat,lims)
    write_xdmf(opt['pfx'],opt['prop'],mat,lims)
    print("Material files generated successfully!")
