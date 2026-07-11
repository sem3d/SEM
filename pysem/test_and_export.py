import os
import numpy as np
import h5py
from gmsh_utils import convert_gmsh_to_h5

def test_integration():
    msh_content = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$PhysicalNames",
        "2",
        '2 20 "bottom_surf"',
        '3 30 "brick_vol"',
        "$EndPhysicalNames",
        "$Nodes",
        "20",
        "1 0.0 0.0 0.0",
        "2 1.0 0.0 0.0",
        "3 1.0 1.0 0.0",
        "4 0.0 1.0 0.0",
        "5 0.0 0.0 1.0",
        "6 1.0 0.0 1.0",
        "7 1.0 1.0 1.0",
        "8 0.0 1.0 1.0",
        "9 0.5 0.0 0.0",
        "10 1.0 0.5 0.0",
        "11 0.5 1.0 0.0",
        "12 0.0 0.5 0.0",
        "13 0.5 0.0 1.0",
        "14 1.0 0.5 1.0",
        "15 0.5 1.0 1.0",
        "16 0.0 0.5 1.0",
        "17 0.0 0.0 0.5",
        "18 1.0 0.0 0.5",
        "19 1.0 1.0 0.5",
        "20 0.0 1.0 0.5",
        "$EndNodes",
        "$Elements",
        "2",
        "1 16 2 20 1 1 2 3 4 9 10 11 12",
        "2 17 2 30 2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20",
        "$EndElements"
    ]
    
    msh_fname = "test_py_mesh.msh"
    h5_fname = "test_py_mesh.h5"
    xmf_fname = "test_py_mesh.xmf"
    
    with open(msh_fname, 'w') as f:
        for line in msh_content:
            f.write(line + "\n")
            
    try:
        print("=== Step 1: Converting GMSH to HDF5/XDMF ===")
        convert_gmsh_to_h5(msh_fname, h5_fname)
        
        # Verify files exist
        assert os.path.exists(h5_fname), "Error: H5 file was not created."
        assert os.path.exists(xmf_fname), "Error: Xdmf file was not created."
        
        print("=== Step 2: Verifying HDF5 file content and shapes ===")
        with h5py.File(h5_fname, 'r') as f_h5:
            # Check Nodes
            assert 'Nodes' in f_h5, "Error: 'Nodes' dataset not found in H5 root."
            nodes_shape = f_h5['/Nodes'].shape
            assert nodes_shape == (28, 3), f"Error: Nodes shape is {nodes_shape}, expected (28, 3)."
            
            # Check /hex27/Hexa27 and /hex27/Mat_Hexa27
            assert 'hex27' in f_h5, "Error: '/hex27' group not found."
            assert 'Hexa27' in f_h5['/hex27'], "Error: 'Hexa27' dataset not found under '/hex27'."
            conn27_shape = f_h5['/hex27/Hexa27'].shape
            assert conn27_shape == (1, 27), f"Error: Hexa27 shape is {conn27_shape}, expected (1, 27)."
            
            mat27_shape = f_h5['/hex27/Mat_Hexa27'].shape
            assert mat27_shape == (1, 1), f"Error: Mat_Hexa27 shape is {mat27_shape}, expected (1, 1)."
            assert f_h5['/hex27/Mat_Hexa27'][0, 0] == 30, "Error: Incorrect material tag for hex27 element."
            
            # Check /surface/Quad9 and /surface/Mat_Quad9
            assert 'surface' in f_h5, "Error: '/surface' group not found."
            assert 'Quad9' in f_h5['/surface'], "Error: 'Quad9' dataset not found under '/surface'."
            conn9_shape = f_h5['/surface/Quad9'].shape
            assert conn9_shape == (1, 9), f"Error: Quad9 shape is {conn9_shape}, expected (1, 9)."
            
            mat9_shape = f_h5['/surface/Mat_Quad9'].shape
            assert mat9_shape == (1, 1), f"Error: Mat_Quad9 shape is {mat9_shape}, expected (1, 1)."
            assert f_h5['/surface/Mat_Quad9'][0, 0] == 20, "Error: Incorrect material tag for quad9 element."
            
            # Check 0-based indexing connectivity
            conn9_data = f_h5['/surface/Quad9'][0]
            assert np.all((conn9_data >= 0) & (conn9_data < 28)), "Error: Quad9 connectivity has out of bounds indices."
            
        print("=== PYTHON INTEGRATION TESTS PASSED SUCCESSFULLY! ===")
        
    finally:
        # Cleanup
        for fname in (msh_fname, h5_fname, xmf_fname):
            if os.path.exists(fname):
                os.remove(fname)

if __name__ == "__main__":
    test_integration()
