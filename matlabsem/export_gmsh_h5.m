function export_gmsh_h5(mesh, h5_filename)
% EXPORT_GMSH_H5  Exporta uma estrutura de malha do Gmsh para arquivos HDF5 e XDMF.
%   Cria h5_filename (.h5) e um arquivo correspondente .xmf no mesmo diretório
%   para leitura direta no ParaView.

[filepath, name, ext] = fileparts(h5_filename);
if isempty(ext)
    h5_filename = fullfile(filepath, [name '.h5']);
    ext = '.h5';
end
xmf_filename = fullfile(filepath, [name '.xmf']);
h5_rel_name = [name ext];

% Se o arquivo H5 já existe, remove para evitar erros de criação de dataset
if exist(h5_filename, 'file')
    delete(h5_filename);
end

% 1. Salvar os nós da malha no root da malha HDF5 (/Nodes)
Nn = size(mesh.nodes, 1);
h5create(h5_filename, '/Nodes', [3, Nn], 'Datatype', 'double');
h5write(h5_filename, '/Nodes', mesh.nodes');

% Tabela de metadados dos tipos de elemento
% { gmsh_type_name, h5_group, h5_dataset, xdmf_type, nnodes }
elem_meta = { ...
    'point',     '/ponto',      'Polyvertex',      'Polyvertex',      1; ...
    'line',      '/linha',      'Seg2',            'Polyline',        2; ...
    'line3',     '/linha',      'Line3',           'Edge_3',          3; ...
    'tri',       '/superficie', 'Tri3',            'Triangle',        3; ...
    'quad',      '/superficie', 'Quad4',           'Quadrilateral',   4; ...
    'tri6',      '/superficie', 'Tri6',            'Triangle_6',      6; ...
    'quad8',     '/superficie', 'Quad8',           'Quadrilateral_8', 8; ...
    'quad9',     '/superficie', 'Quad9',           'Quadrilateral_9', 9; ...
    'tet',       '/volume',     'Tetra4',          'Tetrahedron',     4; ...
    'hex',       '/volume',     'Hexa8',           'Hexahedron',      8; ...
    'tet10',     '/volume',     'Tetra10',         'Tetrahedron_10',  10; ...
    'hex27',     '/hex27',      'Hexa27',          'Hexahedron_27',   27 ...
};

% Agrupar os elementos do mesh.phys por tipo
types = {mesh.phys.type};
unique_types = unique(types);

% Estrutura para guardar informações para o arquivo XDMF
grids_to_write = {};

for it = 1:numel(unique_types)
    type_name = unique_types{it};
    
    % Encontra a linha de metadados correspondente
    meta_row = find(strcmp(elem_meta(:, 1), type_name));
    if isempty(meta_row)
        warning('Tipo de elemento "%s" não reconhecido nos metadados de exportação. Ignorando.', type_name);
        continue;
    end
    
    h5_group = elem_meta{meta_row, 2};
    h5_dataset = elem_meta{meta_row, 3};
    xdmf_type = elem_meta{meta_row, 4};
    nnodes = elem_meta{meta_row, 5};
    
    % Busca elementos e tags de grupos físicos deste tipo
    idx = find(strcmp(types, type_name));
    all_elem = [];
    all_mats = [];
    for i = 1:numel(idx)
        g = idx(i);
        all_elem = [all_elem; mesh.phys(g).elem];
        ne = size(mesh.phys(g).elem, 1);
        all_mats = [all_mats; repmat(mesh.phys(g).tag, ne, 1)];
    end
    
    if isempty(all_elem)
        continue;
    end
    
    Ne = size(all_elem, 1);
    Np = size(all_elem, 2);
    
    if Np ~= nnodes
        error('Erro: Número de nós por elemento (%d) para o tipo %s não condiz com os metadados (%d).', Np, type_name, nnodes);
    end
    
    % Define caminhos de salvamento no HDF5
    h5_path_conn = [h5_group '/' h5_dataset];
    h5_path_mat = [h5_group '/Mat_' h5_dataset];
    
    % Escreve conectividade de elementos (convertendo para 0-based)
    h5create(h5_filename, h5_path_conn, [Np, Ne], 'Datatype', 'int64');
    h5write(h5_filename, h5_path_conn, int64(all_elem' - 1));
    
    % Escreve materiais (tags físicas dos elementos)
    h5create(h5_filename, h5_path_mat, [1, Ne], 'Datatype', 'int64');
    h5write(h5_filename, h5_path_mat, int64(all_mats'));
    
    % Adiciona na lista de grids para escrever no XDMF
    grid.name = type_name;
    grid.xdmf_type = xdmf_type;
    grid.nelems = Ne;
    grid.nnodes = Np;
    grid.path_conn = h5_path_conn;
    grid.path_mat = h5_path_mat;
    grids_to_write{end+1} = grid; %#ok<AGROW>
end

% 2. Escrever o arquivo XDMF (.xmf) correspondente
fid = fopen(xmf_filename, 'w');
if fid < 0
    error('Não foi possível criar o arquivo XDMF: %s', xmf_filename);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '<?xml version="1.0" ?>\n');
fprintf(fid, '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n');
fprintf(fid, '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n');
fprintf(fid, '  <Domain>\n');

% Se tiver múltiplos grids, coloca-os sob uma coleção espacial
if numel(grids_to_write) > 1
    fprintf(fid, '    <Grid Name="mesh" GridType="Collection" CollectionType="Spatial">\n');
    indent = '      ';
else
    indent = '    ';
end

for i = 1:numel(grids_to_write)
    g = grids_to_write{i};
    
    fprintf(fid, '%s<Grid Name="%s" GridType="Uniform">\n', indent, g.name);
    
    % Geometria (coordenadas de nós)
    fprintf(fid, '%s  <Geometry Type="XYZ">\n', indent);
    fprintf(fid, '%s    <DataItem Format="HDF" Dimensions="%d 3" NumberType="Float" Precision="8">%s:/Nodes</DataItem>\n', indent, Nn, h5_rel_name);
    fprintf(fid, '%s  </Geometry>\n', indent);
    
    % Topologia (conectividade dos elementos)
    fprintf(fid, '%s  <Topology Type="%s" NumberOfElements="%d">\n', indent, g.xdmf_type, g.nelems);
    fprintf(fid, '%s    <DataItem Format="HDF" Dimensions="%d %d" NumberType="Int" Precision="8">%s:%s</DataItem>\n', indent, g.nelems, g.nnodes, h5_rel_name, g.path_conn);
    fprintf(fid, '%s  </Topology>\n', indent);
    
    % Atributos das células (campo de materiais)
    fprintf(fid, '%s  <Attribute Name="Mat" Center="Cell" AttributeType="Scalar">\n', indent);
    fprintf(fid, '%s    <DataItem Format="HDF" Dimensions="%d 1" NumberType="Int" Precision="8">%s:%s</DataItem>\n', indent, g.nelems, h5_rel_name, g.path_mat);
    fprintf(fid, '%s  </Attribute>\n', indent);
    
    fprintf(fid, '%s</Grid>\n', indent);
end

if numel(grids_to_write) > 1
    fprintf(fid, '    </Grid>\n');
end

fprintf(fid, '  </Domain>\n');
fprintf(fid, '</Xdmf>\n');

fprintf('Arquivos salvos com sucesso:\n  HDF5: %s\n  XDMF: %s\n', h5_filename, xmf_filename);

end
