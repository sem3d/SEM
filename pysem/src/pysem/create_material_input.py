import h5py
import numpy as np
import argparse

from pysem.compute_pml_length import compute_apow, compute_qmu, compute_qkappa


def generate_material_input(h5_filename, output_filename, has_fluid=False,
                             vp_solid=6300.0, vs_solid=4762.0, rho_solid=2000.0,
                             vp_fluid=1500.0, rho_fluid=1000.0,
                             npow=2, rcdb=-80.0,
                             qmu_coeff=0.1, qp_qs_ratio=2.0,
                             no_attenuation=False):
    """
    Processes a 3D hexahedral mesh from an HDF5 file and generates a material.input
    file with calculated PML properties and correct topological tags.

    Parameters:
    -----------
    h5_filename : str
        Path to the input .h5 mesh file (e.g., 'Argostoli_SRTM3_withBorB_0001_0000.h5')
    output_filename : str
        Path to save the generated 'material.input' file
    has_fluid : bool
        Set to True if material tag 1 is explicitly designated as a fluid (ocean) layer.
    vp_solid, vs_solid, rho_solid : float
        P/S-wave speed [m/s] and density [kg/m3] of the solid (crust/bedrock).
    vp_fluid, rho_fluid : float
        P-wave speed [m/s] and density [kg/m3] of the fluid (ocean); Vs is 0 by definition.
    npow : int
        PML polynomial damping order (material.input "npow" column).
    rcdb : float
        Target PML reflection coefficient, in dB (e.g. -80 for 80 dB attenuation).
    qmu_coeff : float
        Coefficient in the Qmu = qmu_coeff*Vs[m/s] rule of thumb.
    qp_qs_ratio : float
        Assumed Qp/Qs ratio used to derive Qkappa from Qmu via compute_qkappa().
    no_attenuation : bool
        If True, write Qkappa=Qmu=0.0 (SEM3D attenuation off) instead of the
        computed values, matching the previous hardcoded behavior.
    """
    # -------------------------------------------------------------------------
    # 1. User-Customizable Physical Parameters
    # -------------------------------------------------------------------------
    if no_attenuation:
        qkappa_solid = qmu_solid = 0.0
        qkappa_fluid = qmu_fluid = 0.0
    else:
        qmu_solid = compute_qmu(vs_solid, coeff=qmu_coeff)
        qkappa_solid = compute_qkappa(vp_solid, vs_solid, qmu_solid, qp_qs_ratio=qp_qs_ratio)
        # Vs_fluid == 0: Qmu is meaningless for a fluid, Qkappa reduces to Qp = qp_qs_ratio*qmu_solid
        qmu_fluid = 0.0
        qkappa_fluid = qp_qs_ratio * qmu_solid

    vs_fluid = 0.0

    # PML standard properties
    p_order = npow
    a_amplitude = compute_apow(npow, rcdb)

    # -------------------------------------------------------------------------
    # 2. Read Mesh Data from HDF5
    # -------------------------------------------------------------------------
    print(f"Reading mesh datasets from {h5_filename}...")
    with h5py.File(h5_filename, 'r') as f:
        nodes = f['Nodes'][:]
        hexa8 = f['Sem3D/Hexa8'][:]
        mat = f['Sem3D/Mat'][:]
        
    all_tags = np.unique(mat)
    print(f"Detected material tags in mesh: {all_tags}")
    
    # -------------------------------------------------------------------------
    # 3. Detect Core Bounding Box Boundaries
    # -------------------------------------------------------------------------
    core_tags = [0, 1] if has_fluid else [0]
    core_elem_mask = np.isin(mat, core_tags)
    core_nodes_idx = np.unique(hexa8[core_elem_mask])
    core_coords = nodes[core_nodes_idx]
    
    x_min_core = np.min(core_coords[:, 0])
    x_max_core = np.max(core_coords[:, 0])
    y_min_core = np.min(core_coords[:, 1])
    y_max_core = np.max(core_coords[:, 1])
    z_min_core = np.min(core_coords[:, 2])
    
    # -------------------------------------------------------------------------
    # 4. Topological PML Layer Classification (Solid vs Fluid PML)
    # -------------------------------------------------------------------------
    # Group unique nodes owned by each material tag
    tag_nodes = {}
    for tag in all_tags:
        elem_mask = (mat == tag)
        tag_nodes[tag] = set(hexa8[elem_mask].flatten())
        
    pml_type = {}  # Maps PML tag -> 0 (Solid PML) or 1 (Fluid PML)
    
    # Initialize tags directly touching the core domains
    # for tag in all_tags:
    #     if tag in core_tags:
    #         continue
    #     if not tag_nodes[tag].isdisjoint(tag_nodes[0]):
    #         pml_type[tag] = 0
    #     elif has_fluid and not tag_nodes[tag].isdisjoint(tag_nodes[1]):
    #         pml_type[tag] = 1
    # Initialize tags directly touching the core domains
    for tag in all_tags:
        if tag in core_tags:
            continue
            
        # Prioritize checking for Fluid contact first
        if has_fluid and not tag_nodes[tag].isdisjoint(tag_nodes[1]):
            pml_type[tag] = 1  # 1 = Fluid PML ('L')
        # If no fluid contact, check for Solid contact
        elif not tag_nodes[tag].isdisjoint(tag_nodes[0]):
            pml_type[tag] = 0  # 0 = Solid PML ('P')
            
    # Propagate topological association outwards to nested/outer PML blocks
    changed = True
    while changed:
        changed = False
        for tag1 in all_tags:
            if tag1 in core_tags or tag1 in pml_type:
                continue
            for tag2 in list(pml_type.keys()):
                if not tag_nodes[tag1].isdisjoint(tag_nodes[tag2]):
                    pml_type[tag1] = pml_type[tag2]
                    changed = True
                    break

    # -------------------------------------------------------------------------
    # 5. Geometrical Interface & Outward Normal Computation
    # -------------------------------------------------------------------------
    tol = 1e-2  # Small tolerance for spatial alignment checks
    pml_properties = {}
    pml_tags = [tag for tag in all_tags if tag not in core_tags]
    
    for tag in pml_tags:
        elem_mask = (mat == tag)
        nodes_idx = np.unique(hexa8[elem_mask])
        coords_m = nodes[nodes_idx]
        
        xmin, xmax = np.min(coords_m[:, 0]), np.max(coords_m[:, 0])
        ymin, ymax = np.min(coords_m[:, 1]), np.max(coords_m[:, 1])
        zmin, zmax = np.min(coords_m[:, 2]), np.max(coords_m[:, 2])
        
        posX, widthX = 0.0, 0.0
        posY, widthY = 0.0, 0.0
        posZ, widthZ = 0.0, 0.0
        
        # Check X Direction (-x or +x outward normal)
        if xmax <= x_min_core + tol:  # PML is on the -X side (solid is on positive X)
            posX = xmax
            widthX = xmin - xmax      # Changed to negative to absorb on -x
        elif xmin >= x_max_core - tol: # PML is on the +X side (solid is on negative X)
            posX = xmin
            widthX = xmax - xmin      # Keeps positive to absorb on +x
            
        # Check Y Direction (-y or +y outward normal)
        if ymax <= y_min_core + tol:  # PML is on the -Y side (solid is on positive Y)
            posY = ymax
            widthY = ymin - ymax      # Changed to negative to absorb on -y
        elif ymin >= y_max_core - tol: # PML is on the +Y side (solid is on negative Y)
            posY = ymin
            widthY = ymax - ymin      # Keeps positive to absorb on +y
        
        # Check Z Direction (-z outward normal for the bottom PML)
        if zmin < z_min_core + tol:   # Robustly detects bottom PML blocks even with irregular topography
            posZ = zmax               # Interface is at the top of the PML block (boundary with core)
            widthZ = zmin - zmax      # Negative width to cleanly absorb on -z
        else:                         # No PML on the upper +Z surface
            posZ = 0.0
            widthZ = 0.0
        pml_properties[tag] = (posX, widthX, posY, widthY, posZ, widthZ)

    # -------------------------------------------------------------------------
    # 6. Write Formatted material.input File
    # -------------------------------------------------------------------------
    print(f"Writing structured output to {output_filename}...")
    with open(output_filename, 'w') as fout:
        # Part 1: Header and Material definitions
        fout.write(f"{len(all_tags)}\n")
        for tag in sorted(all_tags):
            if tag == 0:
                fout.write(f"S {vp_solid:.6f} {vs_solid:.6f} {rho_solid:.6f} {qkappa_solid:.6f} {qmu_solid:.6f}\n")
            elif has_fluid and tag == 1:
                fout.write(f"F {vp_fluid:.6f} {vs_fluid:.6f} {rho_fluid:.6f} {qkappa_fluid:.6f} {qmu_fluid:.6f}\n")
            else:
                m_type = pml_type.get(tag, 0)
                tag_char = "P" if m_type == 0 else "L"
                vp = vp_solid if m_type == 0 else vp_fluid
                vs = vs_solid if m_type == 0 else vs_fluid
                rho = rho_solid if m_type == 0 else rho_fluid
                qk = qkappa_solid if m_type == 0 else qkappa_fluid
                qm = qmu_solid if m_type == 0 else qmu_fluid
                fout.write(f"{tag_char} {vp:.6f} {vs:.6f} {rho:.6f} {qk:.6f} {qm:.6f}\n")
                
        # Part 2: PML Domain Configuration Mapping
        fout.write("# PML properties\n")
        fout.write("# npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat\n")
        for tag in sorted(pml_tags):
            posX, widthX, posY, widthY, posZ, widthZ = pml_properties[tag]
            m_type = pml_type.get(tag, 0)
            fout.write(f"{p_order} {a_amplitude:.6f} {posX:.6f} {widthX:.6f} {posY:.6f} {widthY:.6f} {posZ:.6f} {widthZ:.6f} {m_type}\n")
            
    print("Process complete! material.input file successfully generated.")

def main():
    parser = argparse.ArgumentParser(description="Generate material.input for seismic simulations.")
    parser.add_argument("-i", "--input", help="Path to the .h5 mesh file")
    parser.add_argument("-o", "--output", default="material.input", help="Output filename")
    parser.add_argument("-f", "--fluid", action="store_true", help="Include fluid domain properties")

    parser.add_argument("--vp-solid", type=float, default=6300.0, help="Solid P-wave speed [m/s]")
    parser.add_argument("--vs-solid", type=float, default=4762.0, help="Solid S-wave speed [m/s]")
    parser.add_argument("--rho-solid", type=float, default=2000.0, help="Solid density [kg/m3]")
    parser.add_argument("--vp-fluid", type=float, default=1500.0, help="Fluid P-wave speed [m/s]")
    parser.add_argument("--rho-fluid", type=float, default=1000.0, help="Fluid density [kg/m3]")

    parser.add_argument("--npow", type=int, default=2, help="PML polynomial damping order")
    parser.add_argument("--rcdb", type=float, default=-80.0,
                         help="Target PML reflection coefficient in dB (e.g. -80 for 80 dB)")
    parser.add_argument("--qmu-coeff", type=float, default=0.1,
                         help="Coefficient in the Qmu = qmu_coeff*Vs[m/s] rule of thumb")
    parser.add_argument("--qp-qs-ratio", type=float, default=2.0,
                         help="Assumed Qp/Qs ratio used to derive Qkappa from Qmu")
    parser.add_argument("--no-attenuation", action="store_true",
                         help="Write Qkappa=Qmu=0.0 instead of computed values")

    args = parser.parse_args()

    generate_material_input(args.input, args.output, args.fluid,
                             vp_solid=args.vp_solid, vs_solid=args.vs_solid, rho_solid=args.rho_solid,
                             vp_fluid=args.vp_fluid, rho_fluid=args.rho_fluid,
                             npow=args.npow, rcdb=args.rcdb,
                             qmu_coeff=args.qmu_coeff, qp_qs_ratio=args.qp_qs_ratio,
                             no_attenuation=args.no_attenuation)
    print(f"Successfully generated {args.output} from {args.input}")

if __name__ == "__main__":
    main()