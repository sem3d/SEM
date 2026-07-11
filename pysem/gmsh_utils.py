import os
import re
import numpy as np
import h5py

def gmsh_npe(etype):
    """Returns the number of nodes per Gmsh element type."""
    npe_map = {
        1: 2,   # 2-node line
        2: 3,   # 3-node triangle
        3: 4,   # 4-node quadrangle
        4: 4,   # 4-node tetrahedron
        5: 8,   # 8-node hexahedron
        6: 6,   # 6-node prism
        7: 5,   # 5-node pyramid
        8: 3,   # 3-node 2nd order line
        9: 6,   # 6-node 2nd order triangle
        10: 9,  # 9-node 2nd order quadrangle
        11: 10, # 10-node 2nd order tetrahedron
        12: 27, # 27-node 2nd order hexahedron
        15: 1,  # 1-node point
        16: 8,  # 8-node 2nd order quadrangle
        17: 20, # 20-node 2nd order hexahedron
        18: 15, # 15-node 2nd order prism
        19: 13, # 13-node 2nd order pyramid
    }
    if etype not in npe_map:
        raise ValueError(f"Gmsh element type {etype} is not supported.")
    return npe_map[etype]

def gmsh_dim(etype):
    """Returns the geometric dimension of a Gmsh element type."""
    if etype == 15:
        return 0
    elif etype in (1, 8):
        return 1
    elif etype in (2, 3, 9, 10, 16):
        return 2
    elif etype in (4, 5, 6, 7, 11, 12, 17, 18, 19):
        return 3
    else:
        return -1

def gmsh_type_name(etype):
    """Returns the descriptive string name of a Gmsh element type."""
    names = {
        1: 'line',
        2: 'tri',
        3: 'quad',
        4: 'tet',
        5: 'hex',
        6: 'prism',
        7: 'pyramid',
        8: 'line3',
        9: 'tri6',
        10: 'quad9',
        11: 'tet10',
        12: 'hex27',
        15: 'point',
        16: 'quad8',
        17: 'hex20',
        18: 'prism15',
        19: 'pyramid13'
    }
    return names.get(etype, 'unknown')

def convert_gmsh_to_h5(msh_filename, h5_filename):
    """
    Reads a Gmsh .msh file (version 2.2 or 4.1), converts 2nd order serendipity elements 
    (quad8/hex20) to Lagrangian elements (quad9/hex27), and exports them to HDF5/XDMF.
    """
    
    print(f"Reading Gmsh mesh file: {msh_filename}...")
    if not os.path.exists(msh_filename):
        raise FileNotFoundError(f"Gmsh file not found: {msh_filename}")
        
    version = 2.2
    phys_dim = []
    phys_tag = []
    phys_name = []
    
    point_ep = []
    curve_ep = []
    surf_ep = []
    vol_ep = []
    
    node_tags = []
    nodes = []
    
    max_type = 20
    elem_conn = [None] * (max_type + 1)
    elem_etag = [None] * (max_type + 1)
    for i in range(max_type + 1):
        elem_conn[i] = []
        elem_etag[i] = []
        
    # ----------------------------------------------------
    # 1. READ AND PARSE THE GMSH MSH FILE
    # ----------------------------------------------------
    print("Parsing sections from Gmsh file...")
    with open(msh_filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            
            if line == '$MeshFormat':
                version_info = f.readline().strip().split()
                version = float(version_info[0])
                if version != 2.2 and version != 4.1:
                    raise ValueError(f"Gmsh MSH format version {version} is not supported. Only 2.2 and 4.1 are accepted.")
                    
            elif line == '$PhysicalNames':
                np_names = int(f.readline().strip())
                for _ in range(np_names):
                    name_line = f.readline().strip()
                    m = re.match(r'^(\d+)\s+(-?\d+)\s+"(.*)"', name_line)
                    if m:
                        phys_dim.append(int(m.group(1)))
                        phys_tag.append(int(m.group(2)))
                        phys_name.append(m.group(3))
                        
            elif line == '$Entities':
                counts = [int(x) for x in f.readline().strip().split()]
                nP, nC, nS, nV = counts[0], counts[1], counts[2], counts[3]
                
                # Points
                for _ in range(nP):
                    parts = f.readline().strip().split()
                    et = int(parts[0])
                    m_phys = int(parts[4])
                    for q in range(m_phys):
                        point_ep.append([et, int(parts[5 + q])])
                        
                # Curves
                for _ in range(nC):
                    parts = f.readline().strip().split()
                    et = int(parts[0])
                    m_phys = int(parts[7])
                    for q in range(m_phys):
                        curve_ep.append([et, int(parts[8 + q])])
                        
                # Surfaces
                for _ in range(nS):
                    parts = f.readline().strip().split()
                    et = int(parts[0])
                    m_phys = int(parts[7])
                    for q in range(m_phys):
                        surf_ep.append([et, int(parts[8 + q])])
                        
                # Volumes
                for _ in range(nV):
                    parts = f.readline().strip().split()
                    et = int(parts[0])
                    m_phys = int(parts[7])
                    for q in range(m_phys):
                        vol_ep.append([et, int(parts[8 + q])])
                        
            elif line == '$Nodes':
                print("Parsing nodes block...")
                if version >= 4.0:
                    header = [int(x) for x in f.readline().strip().split()]
                    nB, nN = header[0], header[1]
                    node_tags = np.zeros(nN, dtype=np.int64)
                    nodes = np.zeros((nN, 3), dtype=np.float64)
                    ptr = 0
                    for _ in range(nB):
                        bh = [int(x) for x in f.readline().strip().split()]
                        edim, param, n = bh[0], bh[2], bh[3]
                        
                        tg = []
                        while len(tg) < n:
                            tg.extend([int(x) for x in f.readline().strip().split()])
                        tg = np.array(tg, dtype=np.int64)
                        
                        coords_read = []
                        n_coords = (3 + edim) * n if param != 0 else 3 * n
                        while len(coords_read) < n_coords:
                            coords_read.extend([float(x) for x in f.readline().strip().split()])
                        coords_read = np.array(coords_read, dtype=np.float64)
                        
                        if param == 0:
                            C = coords_read.reshape(n, 3)
                        else:
                            C = coords_read.reshape(n, 3 + edim)[:, :3]
                            
                        node_tags[ptr:ptr+n] = tg
                        nodes[ptr:ptr+n, :] = C
                        ptr += n
                else:
                    nN = int(f.readline().strip())
                    node_tags = np.zeros(nN, dtype=np.int64)
                    nodes = np.zeros((nN, 3), dtype=np.float64)
                    for i in range(nN):
                        parts = f.readline().strip().split()
                        node_tags[i] = int(parts[0])
                        nodes[i, :] = [float(x) for x in parts[1:4]]
                        
            elif line == '$Elements':
                print("Parsing elements block...")
                if version >= 4.0:
                    header = [int(x) for x in f.readline().strip().split()]
                    nB = header[0]
                    for _ in range(nB):
                        bh = [int(x) for x in f.readline().strip().split()]
                        etag, etype, n = bh[1], bh[2], bh[3]
                        nn = gmsh_npe(etype)
                        
                        if etype >= len(elem_conn):
                            new_size = max(len(elem_conn) * 2, etype + 1)
                            elem_conn.extend([[] for _ in range(new_size - len(elem_conn))])
                            elem_etag.extend([[] for _ in range(new_size - len(elem_etag))])
                            
                        vals = []
                        n_vals = (1 + nn) * n
                        while len(vals) < n_vals:
                            vals.extend([int(x) for x in f.readline().strip().split()])
                        V = np.array(vals, dtype=np.int64).reshape(n, 1 + nn)
                        conn = V[:, 1:]
                        
                        elem_conn[etype].append(conn)
                        elem_etag[etype].append(np.repeat(etag, n))
                else:
                    nE = int(f.readline().strip())
                    for _ in range(nE):
                        parts = [int(x) for x in f.readline().strip().split()]
                        etype = parts[1]
                        ntags = parts[2]
                        ptag = parts[3]
                        conn = parts[3 + ntags:]
                        
                        if etype >= len(elem_conn):
                            new_size = max(len(elem_conn) * 2, etype + 1)
                            elem_conn.extend([[] for _ in range(new_size - len(elem_conn))])
                            elem_etag.extend([[] for _ in range(new_size - len(elem_etag))])
                            
                        elem_conn[etype].append(np.array([conn], dtype=np.int64))
                        elem_etag[etype].append(np.array([ptag], dtype=np.int64))
                        
    for i in range(len(elem_conn)):
        if len(elem_conn[i]) > 0:
            elem_conn[i] = np.vstack(elem_conn[i])
            elem_etag[i] = np.concatenate(elem_etag[i])
        else:
            elem_conn[i] = np.zeros((0, 0), dtype=np.int64)
            elem_etag[i] = np.zeros(0, dtype=np.int64)

    # ----------------------------------------------------
    # 2. REMAP NODES TO CONTIGUOUS 0-BASED INDEXING
    # ----------------------------------------------------
    print("Remapping node indices to contiguous 0-based indexing...")
    nodes = np.array(nodes, dtype=np.float64)
    node_tags = np.array(node_tags, dtype=np.int64)
    
    max_tag = np.max(node_tags) if len(node_tags) > 0 else 0
    id2idx = np.zeros(max_tag + 1, dtype=np.int64)
    id2idx[node_tags] = np.arange(len(node_tags))
    
    for etype in range(len(elem_conn)):
        if elem_conn[etype].size > 0:
            elem_conn[etype] = id2idx[elem_conn[etype]]
            
    # ----------------------------------------------------
    # 3. CONVERT SERENDIPITY ELEMENTS TO LAGRANGIAN
    # ----------------------------------------------------
    print("Converting serendipity elements to Lagrangian...")
    # Quad8 (type 16) -> Quad9 (type 10)
    if len(elem_conn) >= 17 and elem_conn[16].shape[0] > 0:
        print("  Converting quad8 elements to quad9...")
        Ne = elem_conn[16].shape[0]
        corners = elem_conn[16][:, :4]
        new_coords = np.mean(nodes[corners], axis=1)
        
        N_old = len(nodes)
        nodes = np.vstack([nodes, new_coords])
        
        new_ids = np.arange(N_old, N_old + Ne, dtype=np.int64)
        node_tags = np.concatenate([node_tags, new_ids + 1])
        
        conn9 = np.hstack([elem_conn[16], new_ids.reshape(-1, 1)])
        
        if elem_conn[10].shape[0] > 0:
            elem_conn[10] = np.vstack([elem_conn[10], conn9])
            elem_etag[10] = np.concatenate([elem_etag[10], elem_etag[16]])
        else:
            elem_conn[10] = conn9
            elem_etag[10] = elem_etag[16]
            
        elem_conn[16] = np.zeros((0, 8), dtype=np.int64)
        elem_etag[16] = np.zeros(0, dtype=np.int64)
        
    # Hex20 (type 17) -> Hex27 (type 12)
    if len(elem_conn) >= 18 and elem_conn[17].shape[0] > 0:
        print("  Converting hex20 elements to hex27 (performing face deduplication)...")
        Ne = elem_conn[17].shape[0]
        
        local_faces = np.array([
            [0, 1, 2, 3], # face 1
            [4, 5, 6, 7], # face 2
            [0, 1, 5, 4], # face 3
            [1, 2, 6, 5], # face 4
            [2, 3, 7, 6], # face 5
            [3, 0, 4, 7]  # face 6
        ], dtype=np.int64)
        
        all_faces_corners = np.zeros((6 * Ne, 4), dtype=np.int64)
        for f in range(6):
            all_faces_corners[f*Ne : (f+1)*Ne, :] = elem_conn[17][:, local_faces[f]]
            
        sorted_faces = np.sort(all_faces_corners, axis=1)
        unique_sorted, inverse_indices = np.unique(sorted_faces, axis=0, return_inverse=True)
        NumUniqueFaces = unique_sorted.shape[0]
        
        face_coords = np.mean(nodes[unique_sorted], axis=1)
        
        N_old = len(nodes)
        nodes = np.vstack([nodes, face_coords])
        new_face_tags = np.arange(N_old, N_old + NumUniqueFaces, dtype=np.int64)
        node_tags = np.concatenate([node_tags, new_face_tags + 1])
        
        element_face_nodes = np.zeros((Ne, 6), dtype=np.int64)
        for f in range(6):
            element_face_nodes[:, f] = N_old + inverse_indices[f*Ne : (f+1)*Ne]
            
        corners = elem_conn[17][:, :8]
        vol_coords = np.mean(nodes[corners], axis=1)
        
        N_after_faces = len(nodes)
        nodes = np.vstack([nodes, vol_coords])
        new_vol_tags = np.arange(N_after_faces, N_after_faces + Ne, dtype=np.int64)
        node_tags = np.concatenate([node_tags, new_vol_tags + 1])
        element_vol_nodes = new_vol_tags.reshape(-1, 1)
        
        conn27 = np.hstack([elem_conn[17], element_face_nodes, element_vol_nodes])
        
        if elem_conn[12].shape[0] > 0:
            elem_conn[12] = np.vstack([elem_conn[12], conn27])
            elem_etag[12] = np.concatenate([elem_etag[12], elem_etag[17]])
        else:
            elem_conn[12] = conn27
            elem_etag[12] = elem_etag[17]
            
        elem_conn[17] = np.zeros((0, 20), dtype=np.int64)
        elem_etag[17] = np.zeros(0, dtype=np.int64)

    # ----------------------------------------------------
    # 4. PREPARE PHYSICAL GROUP DATA FOR EXPORT
    # ----------------------------------------------------
    print("Grouping physical entities...")
    point_ep = np.array(point_ep, dtype=np.int64) if len(point_ep) > 0 else np.zeros((0, 2), dtype=np.int64)
    curve_ep = np.array(curve_ep, dtype=np.int64) if len(curve_ep) > 0 else np.zeros((0, 2), dtype=np.int64)
    surf_ep  = np.array(surf_ep,  dtype=np.int64) if len(surf_ep)  > 0 else np.zeros((0, 2), dtype=np.int64)
    vol_ep   = np.array(vol_ep,   dtype=np.int64) if len(vol_ep)   > 0 else np.zeros((0, 2), dtype=np.int64)
    
    mesh_phys = []
    
    for ip in range(len(phys_tag)):
        pdim = phys_dim[ip]
        ptag = phys_tag[ip]
        pname = phys_name[ip]
        
        if version >= 4.0:
            if pdim == 0:
                ents = point_ep[point_ep[:, 1] == ptag, 0]
            elif pdim == 1:
                ents = curve_ep[curve_ep[:, 1] == ptag, 0]
            elif pdim == 2:
                ents = surf_ep[surf_ep[:, 1] == ptag, 0]
            elif pdim == 3:
                ents = vol_ep[vol_ep[:, 1] == ptag, 0]
            else:
                ents = []
        else:
            ents = [ptag]
            
        for etype in range(len(elem_conn)):
            if elem_conn[etype].size == 0:
                continue
            if gmsh_dim(etype) != pdim:
                continue
                
            sel = np.isin(elem_etag[etype], ents)
            conn = elem_conn[etype][sel, :]
            if conn.size == 0:
                continue
                
            mesh_phys.append({
                'name': pname,
                'dim': pdim,
                'tag': ptag,
                'type': gmsh_type_name(etype),
                'elem': conn
            })
            
    print(f"Mesh loaded (v{version:.1f}): {nodes.shape[0]} nodes, {len(mesh_phys)} physical groups.")
    for g in mesh_phys:
        print(f"  [{g['type']}] {g['name']:-<12} dim={g['dim']}  {g['elem'].shape[0]:>6} elements")

    # ----------------------------------------------------
    # 5. WRITE TO HDF5 AND XDMF
    # ----------------------------------------------------
    filepath, name_ext = os.path.split(h5_filename)
    name, ext = os.path.splitext(name_ext)
    if not ext:
        h5_filename = os.path.join(filepath, name + '.h5')
        ext = '.h5'
    xmf_filename = os.path.join(filepath, name + '.xmf')
    h5_rel_name = name + ext
    
    print(f"Writing datasets to HDF5 file: {h5_filename}...")
    if os.path.exists(h5_filename):
        os.remove(h5_filename)
        
    elem_meta = {
        'point':     ('/point',      'Polyvertex',      'Polyvertex',      1),
        'line':      ('/line',       'Seg2',            'Polyline',        2),
        'line3':     ('/line',       'Line3',           'Edge_3',          3),
        'tri':       ('/surface',    'Tri3',            'Triangle',        3),
        'quad':      ('/surface',    'Quad4',           'Quadrilateral',   4),
        'tri6':      ('/surface',    'Tri6',            'Triangle_6',      6),
        'quad8':     ('/surface',    'Quad8',           'Quadrilateral_8', 8),
        'quad9':     ('/surface',    'Quad9',           'Quadrilateral_9', 9),
        'tet':       ('/volume',     'Tetra4',          'Tetrahedron',     4),
        'hex':       ('/volume',     'Hexa8',           'Hexahedron',      8),
        'tet10':     ('/volume',     'Tetra10',         'Tetrahedron_10',  10),
        'hex27':     ('/hex27',      'Hexa27',          'Hexahedron_27',   27)
    }
    
    mesh_phys_types = [g['type'] for g in mesh_phys]
    unique_types = set(mesh_phys_types)
    
    grids_to_write = []
    
    with h5py.File(h5_filename, 'w') as f_h5:
        f_h5.create_dataset('/Nodes', data=nodes, dtype='float64')
        
        for type_name in unique_types:
            if type_name not in elem_meta:
                print(f"Warning: Element type '{type_name}' not recognized. Skipping HDF5 write.")
                continue
                
            h5_group, h5_dataset, xdmf_type, nnodes = elem_meta[type_name]
            
            all_elem = []
            all_mats = []
            for g in mesh_phys:
                if g['type'] == type_name:
                    all_elem.append(g['elem'])
                    all_mats.append(np.repeat(g['tag'], g['elem'].shape[0]))
                    
            if len(all_elem) == 0:
                continue
                
            all_elem = np.vstack(all_elem)
            all_mats = np.concatenate(all_mats)
            
            Ne = all_elem.shape[0]
            Np = all_elem.shape[1]
            if Np != nnodes:
                raise ValueError(f"Error: Element connectivity column count ({Np}) for {type_name} does not match nnodes metadata ({nnodes}).")
                
            h5_path_conn = f"{h5_group}/{h5_dataset}"
            h5_path_mat = f"{h5_group}/Mat_{h5_dataset}"
            
            f_h5.create_dataset(h5_path_conn, data=all_elem, dtype='int64')
            f_h5.create_dataset(h5_path_mat, data=all_mats.reshape(-1, 1), dtype='int64')
            
            grids_to_write.append({
                'name': type_name,
                'xdmf_type': xdmf_type,
                'nelems': Ne,
                'nnodes': Np,
                'path_conn': h5_path_conn,
                'path_mat': h5_path_mat
            })
            
    print(f"Writing XDMF file: {xmf_filename}...")
    Nn = nodes.shape[0]
    with open(xmf_filename, 'w') as f_xmf:
        f_xmf.write('<?xml version="1.0" ?>\n')
        f_xmf.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f_xmf.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f_xmf.write('  <Domain>\n')
        
        if len(grids_to_write) > 1:
            f_xmf.write('    <Grid Name="mesh" GridType="Collection" CollectionType="Spatial">\n')
            indent = '      '
        else:
            indent = '    '
            
        for g in grids_to_write:
            f_xmf.write(f'{indent}<Grid Name="{g["name"]}" GridType="Uniform">\n')
            
            f_xmf.write(f'{indent}  <Geometry Type="XYZ">\n')
            f_xmf.write(f'{indent}    <DataItem Format="HDF" Dimensions="{Nn} 3" NumberType="Float" Precision="8">{h5_rel_name}:/Nodes</DataItem>\n')
            f_xmf.write(f'{indent}  </Geometry>\n')
            
            f_xmf.write(f'{indent}  <Topology Type="{g["xdmf_type"]}" NumberOfElements="{g["nelems"]}">\n')
            f_xmf.write(f'{indent}    <DataItem Format="HDF" Dimensions="{g["nelems"]} {g["nnodes"]}" NumberType="Int" Precision="8">{h5_rel_name}:{g["path_conn"]}</DataItem>\n')
            f_xmf.write(f'{indent}  </Topology>\n')
            
            f_xmf.write(f'{indent}  <Attribute Name="Mat" Center="Cell" AttributeType="Scalar">\n')
            f_xmf.write(f'{indent}    <DataItem Format="HDF" Dimensions="{g["nelems"]} 1" NumberType="Int" Precision="8">{h5_rel_name}:{g["path_mat"]}</DataItem>\n')
            f_xmf.write(f'{indent}  </Attribute>\n')
            
            f_xmf.write(f'{indent}</Grid>\n')
            
        if len(grids_to_write) > 1:
            f_xmf.write('    </Grid>\n')
            
        f_xmf.write('  </Domain>\n')
        f_xmf.write('</Xdmf>\n')
        
    print("Files written successfully:")
    print(f"  HDF5: {h5_filename}")
    print(f"  XDMF: {xmf_filename}")
