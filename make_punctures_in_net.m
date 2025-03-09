clear all
clc
tic

d = 1.4 ; % fiber diameter
lx = 750; ly = 750; rho_nodal = 0.02 ;
edge_thkness = 20/ly ;

%% load unfilled network (Change to your directory): 
P = 'C:\Users\mainak\Downloads\DFN_Abaqus\';
F = ['2DNetwork_lf32_d0p02_750x750_Seed_1030'] ;
new = fullfile(P,F) ;

%% the network mat file with nodes and fibers:
load(('C:\Users\mainak\Downloads\DFN_Abaqus\2DNetwork_lf32_d0p02_750x750_Seed_1030\2DNetwork_lf32_d0p02_750x750_Seed_1030.mat')) ;
% no midnodes in fibers at this moment
nodes_set_final = [(1:1:size(nodes, 1)).' nodes] ;
el_set_final = [(1:1:size(fibers, 1)).' fibers] ;
[average_nodal_connectivity] = find_avg_node_connectivity(el_set_final, nodes_set_final, edge_thkness)    % average nodal connectivity, leave edges

[l_fiber] = distance_finder(nodes_set_final, el_set_final, edge_thkness)  % get average fiber length (32.4355 for S3)
SR_circular = d^2/(16*(l_fiber)^2)

% return
%% make the folder that will contain the MAT file of network
date = datestr(now,'mm_dd_yyyy') ;
P = 'C:\Users\mainak\Downloads\DFN_Abaqus\';
F = ['NetSsS5_C4Rh_8EpF_',num2str(lx),'x',num2str(ly),'_rho_',num2str(rho_nodal),'_d_',num2str(d),'_date_', date] ; % in the MATLAB current folder
mkdir(P,F)
new = fullfile(P,F) ;


%% section 7A: punching hole:
[l_fiber] = distance_finder(nodes_set_final, el_set_final, edge_thkness)
fl = l_fiber ;

D = 65 ;
aa = D/2 ;
bb = aa ;

% Specify locations of FOUR circular punctures on the 2D matrix...
x0(1,1) = (( max(nodes_set_final(:, 2)) - min(nodes_set_final(:, 2)) ) / 2) - (lx/2) ;
y0(1,1) = ( max(nodes_set_final(:, 3)) - min(nodes_set_final(:, 3)) ) / 2 - (ly/2) - 55.25*1 ;
x0(2,1) = (( max(nodes_set_final(:, 2)) - min(nodes_set_final(:, 2)) ) / 2) - (lx/2) ;
y0(2,1) = ( max(nodes_set_final(:, 3)) - min(nodes_set_final(:, 3)) ) / 2 - (ly/2) + 55.25*1 ;

% 
x0(3,1) = (( max(nodes_set_final(:, 2)) - min(nodes_set_final(:, 2)) ) / 2) - (lx/2) - 95.7*1 ;
y0(3,1) = ( max(nodes_set_final(:, 3)) - min(nodes_set_final(:, 3)) ) / 2 - (ly/2) ;
x0(4,1) = (( max(nodes_set_final(:, 2)) - min(nodes_set_final(:, 2)) ) / 2) - (lx/2) + 95.7*1 ;
y0(4,1) = ( max(nodes_set_final(:, 3)) - min(nodes_set_final(:, 3)) ) / 2 - (ly/2) ;


for tt = 1:size(x0, 1)
[nodes_set_final, el_set_final] = punch_ellipse_2D_lat_net_rev_VOR(aa, aa, x0(tt), y0(tt), nodes_set_final, el_set_final, fl) ;
end


[l_fiber] = distance_finder(nodes_set_final, el_set_final, edge_thkness)  % get average fiber length 
SR_circular = d^2/(16*(l_fiber)^2)

ly_e = ly * (1-(2*edge_thkness)) ;
[effective_node_nos, rho_nodal_act] = rhonodal_finder(nodes_set_final, el_set_final, edge_thkness, lx, ly_e, 0) ;
rho_nodal_act

%% It adds a midnode in each fiber element 
 [el_set_final, nodes_set_final] = add_midnode_in_fiber(el_set_final, nodes_set_final) ;
% [el_set_final, nodes_set_final] = add_midnode_in_fiber(el_set_final, nodes_set_final) ; % necessary if finer mesh is desired...
% [el_set_final, nodes_set_final] = add_midnode_in_fiber(el_set_final, nodes_set_final) ; 


%% save the fibers and nodes of the network in a MAT file.
Path_lat = char(strcat(new,'\SA2Dnet_data_',num2str(lx),'x',num2str(ly),'_rho_',num2str(rho_nodal),'.mat')) ;
save(Path_lat, 'nodes_set_final' , 'el_set_final', 'edge_thkness', 'l_fiber', 'SR_circular', 'rho_nodal_act', 'D', 'x0', 'y0') ; 

toc

return

