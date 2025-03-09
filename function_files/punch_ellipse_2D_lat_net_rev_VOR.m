% This function file will create an elliptical punch hole at the given location 
% (x0, y0) of the 2D lattice network.
% major, minor radius of the hole is (a, b).
% WRITTEN BY MAINAK SARKAR, UW-MADISON, 2022

function [nodes_set_final, el_set_final] = punch_ellipse_2D_lat_net_rev_VOR(a, b, x0, y0, nodes_set_final, el_set_final, fl)

% el_set_final = [ el_set_final(:, 1) el_set_final(:, 2) el_set_final(:,4) ] ;
fac = -0.01 ; % -0.01
% fac = 0 ; % -0.01

nodes = nodes_set_final ;
mn = size(nodes, 1) ;
nodes = [nodes(:,1:3) , zeros(mn, 1) ] ;
fibers = [el_set_final(:,2), el_set_final(:,3)] ;


% finding all the nodes inside the cylinder
cnt = 0 ;
for i = 1 : size(nodes, 1)
    x = nodes(i, 2) ;
    y = nodes(i, 3) ;
    z = nodes(i, 4) ;
    r_nod = sqrt( (x-x0)^2/a^2 + (y-y0)^2/b^2 ) ;
    if (r_nod-1) < fac*fl
        cnt = cnt + 1 ;
        in_nodes(cnt) = i ;
    end
end

if exist('in_nodes','var')
% remove fibers whose both end nodes are inside the cylinder
art = 0 ;
for j = 1 : size(fibers, 1)
    if ismember(fibers(j,1), in_nodes) == 1 && ismember(fibers(j,2), in_nodes) == 1
    art = art + 1 ;
    nodes_2_delete(art, 1:2) = [fibers(j, 1), fibers(j, 2)] ;
    fibers(j, :) = NaN(1, 2) ; 
    end
end

if exist('nodes_2_delete','var')
nodes_2_delete = nodes_2_delete(:) ;
end 
end

fibers(any(isnan(fibers), 2), :) = [];

if exist('in_nodes','var')
% identify fibers intersecting the cylinder surface and curtail them till
% the surface of the cylinder
cnt2 = 0 ;
for k = 1 : size(fibers, 1)
    n1 = fibers(k, 1) ;
    n2 = fibers(k, 2) ;
    if ismember(n1, in_nodes) == 1 || ismember(n2, in_nodes) == 1
    cnt2 = cnt2 + 1 ;
    xx1 = nodes(n1, 2) ; yy1 = nodes(n1, 3) ; zz1 = nodes(n1, 4) ; 
    xx2 = nodes(n2, 2) ; yy2 = nodes(n2, 3) ; zz2 = nodes(n2, 4) ; 
    mag = sqrt( (xx1 - xx2)^2 + (yy1 - yy2)^2 + (zz1 - zz2)^2 ) ;  
    dir_vec = [ (xx2 - xx1) (yy2 - yy1) (zz2 - zz1) ] ./ mag ;
    nN = 10000 ; 
    tT = linspace(0,1,nN+1)' ;
    intpts = (1-tT) * [xx1 yy1 zz1] + tT * [xx2 yy2 zz2] ;
    for sd = 1 : size(intpts,1)
        rr(sd) = sqrt( (intpts(sd,1)-x0)^2/a^2 + (intpts(sd,2)-y0)^2/b^2 ) ;
    end
        [dsd, index] = min ( ( rr - 1 )) ;
        ttt = rr - 1 ;
        index = knnsearch(ttt',0);
    node_on_surf = intpts(index, :) ;
    if ismember(n1, in_nodes) == 1
    nodes(n1, 2:4) = node_on_surf ;
    surface_node_index(cnt2) = n1 ;
    elseif ismember(n2, in_nodes) == 1
    nodes(n2, 2:4) = node_on_surf ;  
    surface_node_index(cnt2) = n2 ;
    end
    end
end
end

%{
% remove through fibers (fibers passing through the cylinder with no nodes inside cylinder)
for tt = 1 : size(fibers, 1)
    n1 = fibers(tt, 1) ;
    n2 = fibers(tt, 2) ;
%     if (ismember(n1, in_nodes) == 0 && ismember(n2, in_nodes) == 0) 
    xx1 = nodes(n1, 2) ; yy1 = nodes(n1, 3) ; zz1 = nodes(n1, 4) ; 
    xx2 = nodes(n2, 2) ; yy2 = nodes(n2, 3) ; zz2 = nodes(n2, 4) ; 
    nN = 10000 ; 
    tT1 = linspace(0,1,nN+1)' ;
    intpts = (1-tT1) * [xx1 yy1 zz1] + tT1 * [xx2 yy2 zz2] ;
    norm_r = sqrt((intpts(:,1)-x0).^2/a^2 + (intpts(:,2)-y0).^2/b^2) ; 
    for uu = 1 : length(norm_r)
    if norm_r(uu)-1 < fac*fl  
    fibers(tt,:) = NaN(1, 2) ; 
    end
    end
%     end
end
%}

fibers(any(isnan(fibers), 2), :) = [];

% remove nodes inside the cylinder:
if exist('nodes_2_delete','var')

for ij = 1 : length(nodes_2_delete)
    logic_no = (any(surface_node_index(:) == nodes_2_delete(ij))) ;
    if logic_no == 0
    nodes(nodes_2_delete(ij), :) = NaN ;
    end
end

end

nodes(any(isnan(nodes), 2), :) = [];

cnty_mat = [(1:size(nodes,1))' nodes(:,1)] ;

for ab = 1 : size(fibers, 1)
    
ind1 = find(cnty_mat(:, 2) == fibers(ab,1)) ;
ind2 = find(cnty_mat(:, 2) == fibers(ab,2)) ;
    fibers(ab,1) = cnty_mat(ind1,1) ;
    fibers(ab,2) = cnty_mat(ind2,1) ;
    
end

 
nodes = nodes(:,2:4) ;

[fibers] = zero_len_fiber_remover_3D(nodes, fibers) ;
[fibers, nodes] = excess_node_remover_3D(fibers, nodes) ;



el_set_final = [(1:size(fibers,1))', fibers(:,1), fibers(:,2)] ;
nodes_set_final = [(1:size(nodes,1))'  , nodes(:,1:2) ] ;

return


