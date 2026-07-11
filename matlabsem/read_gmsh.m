function mesh = read_gmsh(fname)
% READ_GMSH_PHYS  Lê malha Gmsh (.msh 2.2 ou 4.1) e separa por grupos físicos.
%   Suporta elementos de ordem 1 e ordem 2, com conversão automática de
%   elementos serendipity (quad8 e hex20) para Lagrangianos (quad9 e hex27).
%
%   mesh.nodes        [Nn x 3]  coordenadas dos nós (numeração contígua 1..Nn)
%   mesh.node_ids     [Nn x 1]  tags originais dos nós no Gmsh
%   mesh.phys(g).name           nome do grupo físico
%   mesh.phys(g).dim            dimensão geométrica (0=ponto, 1=curva/linha, 2=superfície, 3=volume)
%   mesh.phys(g).tag            tag físico
%   mesh.phys(g).type           tipo do elemento ('point', 'line', 'line3', 'tri', 'quad',
%                               'tri6', 'quad9', 'tet', 'hex', 'tet10', 'hex27', etc.)
%   mesh.phys(g).elem [Ne x Np] conectividade dos elementos (índices em mesh.nodes)

fid = fopen(fname,'r');
if fid < 0, error('Não foi possível abrir %s', fname); end
cleanup = onCleanup(@() fclose(fid));

version   = 2.2;
phys_dim  = []; phys_tag = []; phys_name = {};
pointEP   = zeros(0,2);   % [entidade_ponto,      tag_fisico]
curveEP   = zeros(0,2);   % [entidade_curva,      tag_fisico]
surfEP    = zeros(0,2);   % [entidade_superficie, tag_fisico]
volEP     = zeros(0,2);   % [entidade_volume,     tag_fisico]
node_tags = []; nodes = [];
max_type  = 20;
elem_conn = cell(max_type, 1);
elem_etag = cell(max_type, 1);

while ~feof(fid)
    line = strtrim(fgetl(fid));
    switch line
        case '$MeshFormat'
            v = sscanf(fgetl(fid),'%f',1);
            version = v;
            if version ~= 2.2 && version ~= 4.1
                error('Versão do formato Gmsh MSH (%.1f) não suportada. Apenas 2.2 e 4.1 são suportadas.', version);
            end

        case '$PhysicalNames'
            np = str2double(fgetl(fid));
            for i = 1:np
                tk = regexp(strtrim(fgetl(fid)), ...
                            '^(\d+)\s+(-?\d+)\s+"(.*)"','tokens','once');
                phys_dim(end+1,1)  = str2double(tk{1}); %#ok<*AGROW>
                phys_tag(end+1,1)  = str2double(tk{2});
                phys_name{end+1,1} = tk{3};
            end

        case '$Entities'   % apenas formato 4.x
            c = sscanf(fgetl(fid),'%d',4);  % nP nC nS nV
            nP=c(1); nC=c(2); nS=c(3); nV=c(4);
            for i=1:nP                       % pontos: tag x y z nPhys phys...
                et=fscanf(fid,'%d',1); fscanf(fid,'%f',3);
                m=fscanf(fid,'%d',1); ph=[]; if m>0, ph=fscanf(fid,'%d',m); end
                for q=1:numel(ph), pointEP(end+1,:)=[et ph(q)]; end
            end
            for i=1:nC                       % curvas
                et=fscanf(fid,'%d',1); fscanf(fid,'%f',6);
                m=fscanf(fid,'%d',1); ph=[]; if m>0, ph=fscanf(fid,'%d',m); end
                nb=fscanf(fid,'%d',1); if nb>0, fscanf(fid,'%d',nb); end
                for q=1:numel(ph), curveEP(end+1,:)=[et ph(q)]; end
            end
            for i=1:nS                       % superfícies
                et=fscanf(fid,'%d',1); fscanf(fid,'%f',6);
                m=fscanf(fid,'%d',1); ph=[]; if m>0, ph=fscanf(fid,'%d',m); end
                nb=fscanf(fid,'%d',1); if nb>0, fscanf(fid,'%d',nb); end
                for q=1:numel(ph), surfEP(end+1,:)=[et ph(q)]; end
            end
            for i=1:nV                       % volumes
                et=fscanf(fid,'%d',1); fscanf(fid,'%f',6);
                m=fscanf(fid,'%d',1); ph=[]; if m>0, ph=fscanf(fid,'%d',m); end
                nb=fscanf(fid,'%d',1); if nb>0, fscanf(fid,'%d',nb); end
                for q=1:numel(ph), volEP(end+1,:)=[et ph(q)]; end
            end
            fgetl(fid);                      % consome resto da linha

        case '$Nodes'
            if version >= 4
                h = sscanf(fgetl(fid),'%d',4); nB=h(1); nN=h(2);
                node_tags = zeros(nN,1); nodes = zeros(nN,3); ptr=0;
                for b=1:nB
                    bh = sscanf(fgetl(fid),'%d',4);
                    edim=bh(1); param=bh(3); n=bh(4);
                    tg = fscanf(fid,'%d',n);
                    if param==0
                        C = reshape(fscanf(fid,'%f',3*n),3,n)';
                    else
                        C = reshape(fscanf(fid,'%f',(3+edim)*n),3+edim,n)';
                        C = C(:,1:3);
                    end
                    node_tags(ptr+1:ptr+n)=tg; nodes(ptr+1:ptr+n,:)=C; ptr=ptr+n;
                    fgetl(fid);
                end
            else  % 2.2
                nN = str2double(fgetl(fid));
                D = fscanf(fid,'%d %f %f %f',[4 nN])';
                node_tags = D(:,1); nodes = D(:,2:4); fgetl(fid);
            end

        case '$Elements'
            if version >= 4
                h = sscanf(fgetl(fid),'%d',4); nB=h(1);
                for b=1:nB
                    bh = sscanf(fgetl(fid),'%d',4);
                    etag=bh(2); etype=bh(3); n=bh(4);
                    nn = gmsh_npe(etype);
                    V = reshape(fscanf(fid,'%d',(1+nn)*n),1+nn,n)';
                    conn = V(:,2:end);
                    if etype > max_type
                        new_max = max(max_type * 2, etype);
                        elem_conn{new_max} = [];
                        elem_etag{new_max} = [];
                        max_type = new_max;
                    end
                    elem_conn{etype} = [elem_conn{etype}; conn];
                    elem_etag{etype} = [elem_etag{etype}; repmat(etag,n,1)];
                    fgetl(fid);
                end
            else  % 2.2: elmTag elmType numTags tags... nodes...
                nE = str2double(fgetl(fid));
                for e=1:nE
                    v = sscanf(fgetl(fid),'%d')';
                    etype=v(2); ntags=v(3);
                    ptag = v(4);                       % 1º tag = físico
                    conn = v(4+ntags:end);
                    if etype > max_type
                        new_max = max(max_type * 2, etype);
                        elem_conn{new_max} = [];
                        elem_etag{new_max} = [];
                        max_type = new_max;
                    end
                    elem_conn{etype} = [elem_conn{etype}; conn];
                    elem_etag{etype} = [elem_etag{etype}; ptag];
                end
            end
    end
end

%% --- Remapeia nós para numeração contígua 1..Nn ---
maxtag = max(node_tags);
id2idx = zeros(maxtag,1);
id2idx(node_tags) = (1:numel(node_tags))';
for etype = 1:numel(elem_conn)
    if ~isempty(elem_conn{etype})
        elem_conn{etype} = reshape(id2idx(elem_conn{etype}), size(elem_conn{etype}));
    end
end

%% --- Conversão de quad8 (16) para quad9 (10) ---
if numel(elem_conn) >= 16 && ~isempty(elem_conn{16})
    Ne = size(elem_conn{16}, 1);
    new_coords = zeros(Ne, 3);
    for e = 1:Ne
        corners = elem_conn{16}(e, 1:4);
        new_coords(e, :) = mean(nodes(corners, :), 1);
    end
    N_old = size(nodes, 1);
    nodes = [nodes; new_coords];
    new_ids = (1:Ne)' + max(node_tags);
    node_tags = [node_tags; new_ids];
    
    conn9 = [elem_conn{16}, (N_old + (1:Ne))'];
    
    if numel(elem_conn) < 10
        elem_conn{10} = [];
        elem_etag{10} = [];
    end
    elem_conn{10} = [elem_conn{10}; conn9];
    elem_etag{10} = [elem_etag{10}; elem_etag{16}];
    elem_conn{16} = [];
    elem_etag{16} = [];
end

%% --- Conversão de hex20 (17) para hex27 (12) ---
if numel(elem_conn) >= 17 && ~isempty(elem_conn{17})
    Ne = size(elem_conn{17}, 1);
    all_faces_corners = zeros(6 * Ne, 4);
    local_faces = [ 1, 2, 3, 4; ...
                    5, 6, 7, 8; ...
                    1, 2, 6, 5; ...
                    2, 3, 7, 6; ...
                    3, 4, 8, 7; ...
                    4, 1, 5, 8 ];
    for f = 1:6
        idx = (f-1)*Ne + (1:Ne);
        all_faces_corners(idx, :) = elem_conn{17}(:, local_faces(f, :));
    end
    
    sorted_faces = sort(all_faces_corners, 2);
    [unique_sorted, ~, ic] = unique(sorted_faces, 'rows');
    NumUniqueFaces = size(unique_sorted, 1);
    
    face_coords = zeros(NumUniqueFaces, 3);
    for i = 1:NumUniqueFaces
        face_coords(i, :) = mean(nodes(unique_sorted(i, :), :), 1);
    end
    
    N_old = size(nodes, 1);
    nodes = [nodes; face_coords];
    new_face_tags = (1:NumUniqueFaces)' + max(node_tags);
    node_tags = [node_tags; new_face_tags];
    
    element_face_nodes = zeros(Ne, 6);
    for f = 1:6
        idx = (f-1)*Ne + (1:Ne);
        element_face_nodes(:, f) = N_old + ic(idx);
    end
    
    vol_coords = zeros(Ne, 3);
    for e = 1:Ne
        corners = elem_conn{17}(e, 1:8);
        vol_coords(e, :) = mean(nodes(corners, :), 1);
    end
    
    N_after_faces = size(nodes, 1);
    nodes = [nodes; vol_coords];
    new_vol_tags = (1:Ne)' + max(node_tags);
    node_tags = [node_tags; new_vol_tags];
    element_vol_nodes = (N_after_faces + (1:Ne))';
    
    conn27 = [elem_conn{17}, element_face_nodes, element_vol_nodes];
    
    if numel(elem_conn) < 12
        elem_conn{12} = [];
        elem_etag{12} = [];
    end
    elem_conn{12} = [elem_conn{12}; conn27];
    elem_etag{12} = [elem_etag{12}; elem_etag{17}];
    elem_conn{17} = [];
    elem_etag{17} = [];
end

mesh.nodes    = nodes;
mesh.node_ids = node_tags;
mesh.phys     = struct('name',{},'dim',{},'tag',{},'type',{},'elem',{});

%% --- Agrupa por grupo físico ---
for ip = 1:numel(phys_tag)
    pdim = phys_dim(ip); ptag = phys_tag(ip);
    
    if version >= 4
        switch pdim
            case 0
                ents = pointEP(pointEP(:,2)==ptag, 1);
            case 1
                ents = curveEP(curveEP(:,2)==ptag, 1);
            case 2
                ents = surfEP(surfEP(:,2)==ptag, 1);
            case 3
                ents = volEP(volEP(:,2)==ptag, 1);
            otherwise
                ents = [];
        end
    else
        ents = ptag;
    end
    
    for etype = 1:numel(elem_conn)
        if isempty(elem_conn{etype}), continue; end
        if gmsh_dim(etype) ~= pdim, continue; end
        
        sel = ismember(elem_etag{etype}, ents);
        conn = elem_conn{etype}(sel, :);
        if isempty(conn), continue; end
        
        g = numel(mesh.phys)+1;
        mesh.phys(g).name = phys_name{ip};
        mesh.phys(g).dim  = pdim;
        mesh.phys(g).tag  = ptag;
        mesh.phys(g).type = gmsh_type_name(etype);
        mesh.phys(g).elem = conn;
    end
end

fprintf('Malha lida (v%.1f): %d nós, %d grupos físicos.\n', ...
        version, size(nodes,1), numel(mesh.phys));
for g = 1:numel(mesh.phys)
    fprintf('  [%s] %-12s dim=%d  %6d elementos %s\n', ...
            mesh.phys(g).type, mesh.phys(g).name, mesh.phys(g).dim, ...
            size(mesh.phys(g).elem,1), mesh.phys(g).type);
end
end

%% =========================================================
function nn = gmsh_npe(t)   % nº de nós por tipo de elemento Gmsh
switch t
    case 1,  nn=2;   case 2,  nn=3;   case 3,  nn=4;   case 4,  nn=4;
    case 5,  nn=8;   case 6,  nn=6;   case 7,  nn=5;   case 8,  nn=3;
    case 9,  nn=6;   case 10, nn=9;   case 11, nn=10;  case 12, nn=27;
    case 15, nn=1;   case 16, nn=8;   case 17, nn=20;  case 18, nn=15;
    case 19, nn=13;
    otherwise, error('Tipo de elemento Gmsh %d não suportado.', t);
end
end

%% =========================================================
function d = gmsh_dim(t)   % dimensão geométrica do elemento Gmsh
switch t
    case 15, d = 0; % ponto
    case {1, 8}, d = 1; % linha
    case {2, 3, 9, 10, 16}, d = 2; % superfície
    case {4, 5, 6, 7, 11, 12, 17, 18, 19}, d = 3; % volume
    otherwise, d = -1;
end
end

%% =========================================================
function name = gmsh_type_name(t)   % string descritiva do tipo de elemento Gmsh
switch t
    case 1, name = 'line';
    case 2, name = 'tri';
    case 3, name = 'quad';
    case 4, name = 'tet';
    case 5, name = 'hex';
    case 6, name = 'prism';
    case 7, name = 'pyramid';
    case 8, name = 'line3';
    case 9, name = 'tri6';
    case 10, name = 'quad9';
    case 11, name = 'tet10';
    case 12, name = 'hex27';
    case 15, name = 'point';
    case 16, name = 'quad8';
    case 17, name = 'hex20';
    case 18, name = 'prism15';
    case 19, name = 'pyramid13';
    otherwise, name = 'unknown';
end
end