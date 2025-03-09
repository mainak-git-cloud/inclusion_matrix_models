% clc
% clear all
% close all
% WRITTEN BY MAINAK SARKAR, UW MADISON, 2022

function [final_nodes, incl_elements, ex_nodes, extra_fibers] = ellp_incl_gen_upgraded_rev(A, x0, y0, a, b, nodes_set_final, el_set_final)

% get the circumferential nodes on matrix:
for ttj = 1 : size(x0, 1)
count = 0 ;
for ii = 1 : size(nodes_set_final,1)
x_node = nodes_set_final(ii,2) ;
y_node = nodes_set_final(ii,3) ;
r_node = ( ((x_node-x0(ttj))*cos(A(ttj))+(y_node-y0(ttj))*sin(A(ttj))).^2/a^2 + ((x_node-x0(ttj))*sin(A(ttj))-(y_node-y0(ttj))*cos(A(ttj))).^2/b^2 ) ;
if floor(r_node-1) == 0 || ceil(r_node-1) == 0 
% if (r_node-1) <= fac^fl
count = count + 1 ;
index_cylsurf_nodes{ttj}(count) = ii ;
cylsurf_nodes{ttj}(count, :) = nodes_set_final(ii,2:3) ;
end
end
end

rst = size(nodes_set_final, 1) ;
for ttj = 1 : size(x0, 1)
t = pi/16:pi/16:2*pi ;
cnt = 0 ;
for q = 1 : -.1 : .1
cnt = cnt + 1 ;
x{cnt} = q*a * cos(t) ;
y{cnt} = q*b * sin(t) ;
mm = size(x{cnt}, 2) ;
points_for_tri{cnt} = [x{cnt}' y{cnt}'] ;
end

points_for_tri_all = vertcat(points_for_tri{:}) ;

% rotate the nodes within the ellipse:
R_mat = [cos(A(ttj)), -sin(A(ttj)) ; sin(A(ttj)), cos(A(ttj))] ; % rotn by angle A radian counterclockwise
aligned_oriented_patch_nodes = (R_mat * points_for_tri_all.').' ; % oriented_patch_nodes is n x 2 matrix
points_for_tri_all(:, 1) = aligned_oriented_patch_nodes(:, 1) + x0(ttj) ; 
points_for_tri_all(:, 2) = aligned_oriented_patch_nodes(:, 2) + y0(ttj) ; 

T = delaunayTriangulation(points_for_tri_all) ;
nodes_incl{ttj} = [T.Points] ;
ex_nodes{ttj} = nodes_incl{ttj}(1:end,:) ; 
rst = rst + size(ex_nodes{ttj}, 1) ;
ele_incl{ttj} = [T.ConnectivityList] ;
scaling_mat{ttj}(:, 1) = [(1:size(nodes_incl{ttj}, 1))'] ;
scaling_mat{ttj}(:, 2) = [((rst - size(ex_nodes{ttj}, 1) + 1):(rst))'] ;

for uv = 1 : size(ele_incl{ttj}, 1)
    revised_ele_incl{ttj}(uv, 1) = scaling_mat{ttj}(ele_incl{ttj}(uv, 1),2) ;
    revised_ele_incl{ttj}(uv, 2) = scaling_mat{ttj}(ele_incl{ttj}(uv, 2),2) ;
    revised_ele_incl{ttj}(uv, 3) = scaling_mat{ttj}(ele_incl{ttj}(uv, 3),2) ;
end

triplot(T)
hold on
end

ele_inclusion = vertcat(revised_ele_incl{:}) ;
extra_nodes = vertcat(ex_nodes{:}) ;

final_nodes = [nodes_set_final(:, 2:3) ; extra_nodes ] ;

final_nodes = [(1:size(final_nodes, 1))', final_nodes] ;

qr = size(el_set_final, 1) ;
incl_elements = [qr+(1:size(ele_inclusion, 1))', ele_inclusion] ;

past_ele_nos = incl_elements(end, 1) ;

% size(index_cylsurf_nodes)
% size(x0)
% get extra_connectors
for yj = 1 : size(x0, 1)
for yi = 1 : length(index_cylsurf_nodes{yj})
x_dwn = final_nodes(index_cylsurf_nodes{yj}(yi), 2) ;
y_dwn = final_nodes(index_cylsurf_nodes{yj}(yi), 3) ;
idx_dwn = index_cylsurf_nodes{yj}(yi) ;
cc = 0 ;
for pq = 1 : size(nodes_incl{yj}, 1)
x_up = nodes_incl{yj}(pq, 1) ;
y_up = nodes_incl{yj}(pq, 2) ;
idx_up = scaling_mat{yj}(pq, 2) ;
cc = cc + 1 ;
r(cc) = sqrt((x_up-x_dwn)^2 + (y_up-y_dwn)^2) ;
end
[rmin, ridmin] = min(r) ;
idx_up_taken = scaling_mat{yj}(ridmin, 2) ;
past_ele_nos = past_ele_nos + 1 ;
extra_connectors{yj}(yi, :) = [past_ele_nos idx_dwn idx_up_taken] ;
end
end

extra_fibers = vertcat(extra_connectors{:}) ;

return


