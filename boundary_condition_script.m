% 
clear all
clc
d = 1.4 ; % fiber diameter

%% load network/matrix: Change based on saving directory
P = 'C:\Users\mainak\Downloads\DFN_Abaqus\';
F = ['NetSsS5_C4Rh_750x750_rho_0.02_d_1.4_date_08_01_2024'] ;
new = fullfile(P,F) ;

%% This is the mat file containing nodes and fibers of the matrix
load(('C:\Users\mainak\Downloads\DFN_Abaqus\NetSsS5_C4Rh_750x750_rho_0.02_d_1.4_date_08_01_2024\SA2Dnet_data_750x750_rho_0.02.mat')) ;

% return
[l_fiber] = distance_finder(nodes_set_final, el_set_final, edge_thkness) ;  % get average fiber length 
SR_circular = d^2/(16*(l_fiber)^2) ;

kappa = SR_circular/4 % ratio EI/EAL^2... is in O(1e-4).


% return
%%
for radial_strain = [.15] % chose the applied radial strain (contractile) on inclusions; this example shows 15% radial contractile strain 
% for radial_strain = [0.5 0.7 0.9]
DP_strain = 0 ; % shear strain in step 2, keeping output structure of step 1 as the reference
% Define file names of INP and ODB files
rad_strain_str = num2str(100*radial_strain);
DPstrain_str = num2str(100*DP_strain);
file = [new, '\NetSsS5_C4Rh_exHDxgit7_750x750_R2c2_',rad_strain_str,'_NaNDP_strain_',DPstrain_str,'.inp'] ;
F2 = ['NetSsS5_C4Rh_exHDxgit7_750x750_R2c2_pre_rad_',rad_strain_str,'_NaNDP_POST_',DPstrain_str] ;
odb_file_name = F2 ; % This file will be stored in the current folder of MATLAB after creating in Abaqus CAE post-simulation

% inp generation: 
% partial contraction of few incl:
fraction = 1 ; % set b/w 0 & 1
% dipole_thermal_strain = DP_strain ; % dipole
% Lf = 32 ;
% dipole_length = 2*Lf ;
steps = 1 ;
cycle = 1 ;
% get inp file
[final_fibers, final_nodes] = inp_gen_incl_ctr_plus_recv(steps, cycle, fraction, file, d, 100*radial_strain, el_set_final, nodes_set_final, edge_thkness, D, x0, y0) ;


end
return
run_inp_abaqus_odb_creator(file, odb_file_name)


 