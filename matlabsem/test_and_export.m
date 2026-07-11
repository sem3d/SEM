function test_and_export()
    % Add current directory to path
    addpath(pwd);

    % Define mock MSH 2.2 content with 2nd order serendipity elements (quad8, hex20)
    msh22_content = {
        '$MeshFormat'
        '2.2 0 8'
        '$EndMeshFormat'
        '$PhysicalNames'
        '2'
        '2 20 "bottom_surf"'
        '3 30 "brick_vol"'
        '$EndPhysicalNames'
        '$Nodes'
        '20'
        '1 0.0 0.0 0.0'
        '2 1.0 0.0 0.0'
        '3 1.0 1.0 0.0'
        '4 0.0 1.0 0.0'
        '5 0.0 0.0 1.0'
        '6 1.0 0.0 1.0'
        '7 1.0 1.0 1.0'
        '8 0.0 1.0 1.0'
        '9 0.5 0.0 0.0'
        '10 1.0 0.5 0.0'
        '11 0.5 1.0 0.0'
        '12 0.0 0.5 0.0'
        '13 0.5 0.0 1.0'
        '14 1.0 0.5 1.0'
        '15 0.5 1.0 1.0'
        '16 0.0 0.5 1.0'
        '17 0.0 0.0 0.5'
        '18 1.0 0.0 0.5'
        '19 1.0 1.0 0.5'
        '20 0.0 1.0 0.5'
        '$EndNodes'
        '$Elements'
        '2'
        '1 16 2 20 1 1 2 3 4 9 10 11 12'
        '2 17 2 30 2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'
        '$EndElements'
    };

    % Write temp MSH 2.2 file
    msh_fname = 'test_integration.msh';
    fid = fopen(msh_fname, 'w');
    for i = 1:numel(msh22_content)
        fprintf(fid, '%s\n', msh22_content{i});
    end
    fclose(fid);
    
    % Step 1: Read the mesh (which performs conversion to quad9 and hex27)
    fprintf('=== Passo 1: Lendo malha Gmsh ===\n');
    mesh = read_gmsh(msh_fname);
    
    % Step 2: Export to HDF5 & XDMF
    h5_fname = 'test_integration.h5';
    xmf_fname = 'test_integration.xmf';
    fprintf('=== Passo 2: Exportando para HDF5/XDMF ===\n');
    export_gmsh_h5(mesh, h5_fname);
    
    % Step 3: Verify output files exist
    assert(exist(h5_fname, 'file') == 2, 'Erro: Arquivo H5 não foi criado.');
    assert(exist(xmf_fname, 'file') == 2, 'Erro: Arquivo XDMF não foi criado.');
    
    % Step 4: Verify HDF5 file internal structure
    fprintf('=== Passo 3: Verificando estrutura do HDF5 ===\n');
    info = h5info(h5_fname);
    
    % Check for '/Nodes' dataset
    node_dataset_idx = find(strcmp({info.Datasets.Name}, 'Nodes'));
    assert(~isempty(node_dataset_idx), 'Erro: Dataset /Nodes não encontrado no root do H5.');
    node_size = info.Datasets(node_dataset_idx).Dataspace.Size;
    assert(isequal(node_size, [3, 28]), 'Erro: Tamanho do dataset /Nodes incorreto (esperado [3, 28]).');
    
    % Check for '/hex27' group and its datasets
    hex27_grp_idx = find(strcmp({info.Groups.Name}, '/hex27'));
    assert(~isempty(hex27_grp_idx), 'Erro: Grupo /hex27 não encontrado no H5.');
    hex27_grp = info.Groups(hex27_grp_idx);
    
    hex27_conn_idx = find(strcmp({hex27_grp.Datasets.Name}, 'Hexa27'));
    assert(~isempty(hex27_conn_idx), 'Erro: Dataset Hexa27 não encontrado em /hex27.');
    hex27_conn_size = hex27_grp.Datasets(hex27_conn_idx).Dataspace.Size;
    assert(isequal(hex27_conn_size, [27, 1]), 'Erro: Tamanho do dataset /hex27/Hexa27 incorreto (esperado [27, 1]).');
    
    hex27_mat_idx = find(strcmp({hex27_grp.Datasets.Name}, 'Mat_Hexa27'));
    assert(~isempty(hex27_mat_idx), 'Erro: Dataset Mat_Hexa27 não encontrado em /hex27.');
    hex27_mat_size = hex27_grp.Datasets(hex27_mat_idx).Dataspace.Size;
    assert(isequal(hex27_mat_size, [1, 1]), 'Erro: Tamanho do dataset /hex27/Mat_Hexa27 incorreto (esperado [1, 1]).');
    
    % Check for '/superficie' group and its datasets
    surf_grp_idx = find(strcmp({info.Groups.Name}, '/superficie'));
    assert(~isempty(surf_grp_idx), 'Erro: Grupo /superficie não encontrado no H5.');
    surf_grp = info.Groups(surf_grp_idx);
    
    quad9_conn_idx = find(strcmp({surf_grp.Datasets.Name}, 'Quad9'));
    assert(~isempty(quad9_conn_idx), 'Erro: Dataset Quad9 não encontrado em /superficie.');
    quad9_conn_size = surf_grp.Datasets(quad9_conn_idx).Dataspace.Size;
    assert(isequal(quad9_conn_size, [9, 1]), 'Erro: Tamanho do dataset /superficie/Quad9 incorreto (esperado [9, 1]).');
    
    quad9_mat_idx = find(strcmp({surf_grp.Datasets.Name}, 'Mat_Quad9'));
    assert(~isempty(quad9_mat_idx), 'Erro: Dataset Mat_Quad9 não encontrado em /superficie.');
    
    % Cleanup temp files
    delete(msh_fname);
    delete(h5_fname);
    delete(xmf_fname);
    
    fprintf('=== TESTE DE INTEGRAÇÃO H5/XDMF PASSOU COM SUCESSO! ===\n');
end
