import os
import xml.etree.ElementTree as ET

ppath = "/scratch/ldeabreucorrea/CaseMarie/SEMs/Ref"
res_dir = ppath + "/res"
file_path_results = ppath + "/res/results.xmf"

def read_xmf_with_includes(xdmf_file):
    tree = ET.parse(xdmf_file)
    root = tree.getroot()
    mesh_files = []
    for include in root.findall(".//{http://www.w3.org/2001/XInclude}include"):
        href = include.get("href")
        if href:
            mesh_files.append(href)
    return mesh_files

def parse_mesh_file(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    domain = root.find("Domain")
    grid_coll = domain.find("Grid")
    
    top_data_items = [ET.tostring(item, encoding="unicode") for item in grid_coll.findall("DataItem")]
    sub_grids = [ET.tostring(g, encoding="unicode") for g in grid_coll.findall("Grid")]
    
    return top_data_items, sub_grids

def write_filtered_mesh(output_path, top_data_items, sub_grid_str, space_id):
    with open(output_path, "w", encoding="utf-8") as file:
        file.write('<?xml version="1.0" ?>\n')
        file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        file.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        file.write('<Domain>\n')
        file.write(f'<Grid CollectionType="Temporal" GridType="Collection" Name="space.{space_id:04d}">\n')
        
        for item_str in top_data_items:
            file.write(item_str)
            if not item_str.endswith('\n'):
                file.write('\n')
                
        file.write(sub_grid_str)
        if not sub_grid_str.endswith('\n'):
            file.write('\n')
            
        file.write('</Grid>\n')
        file.write('</Domain>\n')
        file.write('</Xdmf>\n')

def write_spatial_collection(output_path, href_list):
    with open(output_path, "w", encoding="utf-8") as file:
        file.write('<?xml version="1.0" ?>\n')
        file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        file.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        file.write('<Domain>\n')
        file.write('<Grid CollectionType="Spatial" GridType="Collection">\n')

        for href in href_list:
            file.write(f'<xi:include href="{href}" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n')

        file.write('</Grid>\n')
        file.write('</Domain>\n')
        file.write('</Xdmf>\n')

def main():
    print(f"Reading master results file: {file_path_results}")
    mesh_files = read_xmf_with_includes(file_path_results)
    print(f"Found {len(mesh_files)} mesh partition files.")

    if not mesh_files:
        print("No mesh files found in results.xmf!")
        return

    # Cache pre-parsed data for all mesh files to optimize performance
    mesh_data = {}
    print("Pre-parsing mesh partition files...")
    for imesh in mesh_files:
        full_path = os.path.join(res_dir, imesh)
        mesh_data[imesh] = parse_mesh_file(full_path)

    # Determine total timesteps from first mesh file
    sample_mesh = mesh_files[0]
    _, sample_sub_grids = mesh_data[sample_mesh]
    num_timesteps = len(sample_sub_grids)
    print(f"Detected {num_timesteps} timesteps per mesh file.")

    for ts in range(1, num_timesteps + 1):
        timestep_idx = ts - 1  # 0-indexed for sub_grids list
        results_out_path = os.path.join(res_dir, f"results_{ts}.xmf")
        print(f"Processing timestep {ts}/{num_timesteps} -> results_{ts}.xmf")
        
        href_list = []
        space_id = 0
        
        for imesh in mesh_files:
            top_data_items, sub_grids = mesh_data[imesh]
            if timestep_idx < len(sub_grids):
                filtered_mesh_name = f"Filtered_{ts}_{imesh}"
                filtered_mesh_path = os.path.join(res_dir, filtered_mesh_name)
                
                write_filtered_mesh(filtered_mesh_path, top_data_items, sub_grids[timestep_idx], space_id)
                href_list.append(filtered_mesh_name)
                space_id += 1
        
        write_spatial_collection(results_out_path, href_list)

    print("Successfully generated all individual timestep XMF files!")

if __name__ == "__main__":
    main()

