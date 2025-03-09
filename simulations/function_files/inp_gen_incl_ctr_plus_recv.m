% define material and generate an Abaqus input file for fibrous matrix
% Written by Mainak Sarkar, University of Illinois Urbana-Champaign, 2024


function [final_fibers, final_nodes] = inp_gen_incl_ctr_plus_recv(NN, cycle, fraction, file, d, radial_strain, final_fibers, final_nodes, edge_thkness, D, x0, y0)

side_length_y = max(final_nodes(:,3)) - min(final_nodes(:,3)) ; 
side_length_x = max(final_nodes(:,2)) - min(final_nodes(:,2)) ; 

% MATERIAL PROPERTIES
E = 1e-6 ; % Set Young's modulus. Ensure wave speed, sqrt(E/density), is close to 1
% and preferebly below 10 to avoid large spread of eigen values and ILL
% CONDITIONING. If ill conditioning occurs, very small time step will be
% necessary and the solver will fail.
nu = 0;   % Poisson's ratio. If we don't give Abaqus a value, it uses 0.
% Fiber cross sectional area
Area_a = (pi/4) * d^2 ;

MOI_a = (Area_a^2)/(4*pi) ;
MOI_b = (Area_a^2)/(4*pi);
J_a = MOI_a + MOI_b ;


dt=8e-2;    % time increment
alpha_r = 0.02 ; % damping
% Density for dynamic analysis; Ensure that the wave speed, sqrt(E/density), is close to 1
% and preferebly below 10 to avoid large spread of eigen values and ILL
% CONDITIONING. If ill conditioning occurs, very small time step will be
% necessary and the solver will fail. Having ratio E/density slightly above 1 is good for a quasi-static system.
density = 3.6e-07 ; % mass density 
duration = 2000;



d = edge_thkness * side_length_y ;
%%============== Finding Edges (do edits if necessary)============================================
edge_x_pos = find(final_nodes(:,2) > (side_length_x/2 - d) ) ;
edge_x_neg = find(final_nodes(:,2) < (-side_length_x/2 + d) ) ;
edge_y_pos = find(final_nodes(:,3) > (side_length_y/2) - d ) ;
edge_y_neg = find(final_nodes(:,3) < -(side_length_y/2) + d) ;


pp = size(final_nodes, 1) ;
%% modification due to inclusion: Ellipses / circles
angles = zeros(size(x0)) ;
[final_nodes, ele_incl, incl_nodes, extra_fibers] = ellp_incl_gen_upgraded_rev(angles, x0, y0, D/2, (D/2), final_nodes, final_fibers) ;


%%
radial_delta = D * radial_strain ;
counter = 0 ;
for ij = 1 : size(x0, 1)
    nodes_particle = incl_nodes{ij} ;
    ppx = size(nodes_particle, 1) ;
    init = pp ;
    pp = pp + ppx ;
    counter2 = 0 ;
    for jk = (init+1) : (pp)
        counter = counter + 1 ;
        part_nodes_idx(counter, 1) = final_nodes(jk, 1) ;
        counter2 = counter2 + 1 ;
        part_nodes_idx2(counter2, 1) = final_nodes(jk, 1) ;
        xxx = final_nodes(jk, 2) ;
        yyy = final_nodes(jk, 3) ;
        x_comp = x0(ij) - xxx ;
        y_comp = y0(ij) - yyy ;
        norm = sqrt(x_comp^2 + y_comp^2) ;
        part_nodes_DISPx(counter, 1) = norm * radial_strain * (x_comp/norm) ;
        part_nodes_DISPy(counter, 1) = norm * radial_strain * (y_comp/norm) ;
    end
    part_nodes_idx2_incl{ij} = part_nodes_idx2 ;
    clear part_nodes_idx2
end

el_per_incl = counter/size(x0, 1) ;

%%
midpoint_nodes = zeros(length(final_fibers),3);

for j=1:length(final_fibers)
    midpt = .5*(final_nodes(final_fibers(j,2),2:3)+final_nodes(final_fibers(j,3),2:3));
    midpoint_nodes(j,:) = [length(final_nodes)+j midpt];
end

% Redo the elements to include the midpoint nodes
final_fibers = [final_fibers(:,1) final_fibers(:,2) midpoint_nodes(:,1) final_fibers(:,3)] ;
final_nodes = [final_nodes; midpoint_nodes];



%% --- WRITE OUTPUT TO DATA FILE ---

fid = fopen(file,'w');

fprintf(fid,'*HEADING\n');
fprintf(fid,'**\n');
fprintf(fid,'** Quasi-2D Inclusion-Matrix Model Definition\n');
fprintf(fid,'**\n');


% --- NODES ---
fprintf(fid,'*NODE, NSET=Nodes\n');
% fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',final_nodes');
fprintf(fid,'%8.0f,%20.12E,%20.12E\n',final_nodes');

fprintf(fid,'*NSET, NSET=XPOS\n');
fprintf(fid,'%8.0f\n',edge_x_pos');
fprintf(fid,'*NSET, NSET=XNEG\n');
fprintf(fid,'%8.0f\n',edge_x_neg');
fprintf(fid,'*NSET, NSET=YPOS\n');
fprintf(fid,'%8.0f\n',edge_y_pos');
fprintf(fid,'*NSET, NSET=YNEG\n');
fprintf(fid,'%8.0f\n',edge_y_neg');

for ik = 1 : size(part_nodes_idx, 1)
fprintf(fid,['*NSET, NSET=CNODES_', num2str(ik),'\n']);
fprintf(fid,'%8.0f\n',part_nodes_idx(ik, 1)) ;
end


% --- ELEMENTS ---
fprintf(fid,'**\n');
fprintf(fid,'** Elements\n');
fprintf(fid,'**\n');

%linear beam is B31 in 3D, B21 in 2D
%Quadratic beam is B32 in 3D, B22 in 2D
fprintf(fid,'*ELEMENT, TYPE=B32, ELSET=EL_fiber_a\n'); % For beam element set a
% fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');

fprintf(fid,'*ELEMENT, TYPE=B31, ELSET=cnt_fibers\n'); % For connector beam elements set
fprintf(fid,'%8.0f,%8.0f,%8.0f\n',extra_fibers');

fprintf(fid,['*ELEMENT, TYPE=CPS3, ELSET=incl_elements\n']); % For 3 noded triangular element (plane stress)
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',ele_incl');


% % --- SECTION ASSIGNMENT AND MATERIAL PROPS ---

fprintf(fid,'**\n');
fprintf(fid,'** Section Assignment and Material Properties\n');
fprintf(fid,'**\n');

% section for element set a
fprintf(fid,'*BEAM GENERAL SECTION, ELSET=EL_fiber_a, DENSITY=%20.10E, SECTION=GENERAL\n',density); 
% No material given (it's defined here). Density is a required parameter.
% Poisson's ratio is optional; leaving blank sets it to zero.
fprintf(fid,'%20.10E,%20.10E,%20.10E,%20.10E,%20.10E\n',Area_a,MOI_a,0,MOI_b,J_a); % Area | I11 | I12 | I22 | J
fprintf(fid,'\n'); % blank line to accept default direction consines
fprintf(fid,'%20.10E,%20.10E\n',E,E/(2*(1+nu))); % Young's modulus | shear modulus | other entries unneeded
fprintf(fid,'*SECTION POINTS\n'); % Section points are needed to compute stress and strain within element
fprintf(fid,'%10.4f,%10.4f\n',0.5,0); % This gives x1 and x2 position of section points
fprintf(fid,'*DAMPING, ALPHA=%20.10E\n',alpha_r);
% end of section assignment for element set a


% section for connector elements
fprintf(fid,'*BEAM GENERAL SECTION, ELSET=cnt_fibers, DENSITY=%20.10E, SECTION=GENERAL\n',density); 
% No material given (it's defined here). Density is a required parameter.
% Poisson's ratio is optional; leaving blank sets it to zero.
fprintf(fid,'%20.10E,%20.10E,%20.10E,%20.10E,%20.10E\n',Area_a,MOI_a,0,MOI_b,J_a); % Area | I11 | I12 | I22 | J
fprintf(fid,'\n'); % blank line to accept default direction consines
fprintf(fid,'%20.10E,%20.10E\n',10*E,(10*E)/(2*(1+nu))); % Young's modulus | shear modulus | other entries unneeded
fprintf(fid,'*SECTION POINTS\n'); % Section points are needed to compute stress and strain within element
fprintf(fid,'%10.4f,%10.4f\n',0.5,0); % This gives x1 and x2 position of section points
fprintf(fid,'*DAMPING, ALPHA=%20.10E\n',alpha_r);
% end of section assignment for connector elements


fprintf(fid,'** === === Inclusion Material Properties === ===\n');
fprintf(fid,'*MATERIAL, NAME=ThermalParticle\n');
fprintf(fid,'*ELASTIC, TYPE=ISOTROPIC\n');
fprintf(fid,'%20.10E,%20.10E,%20.10E\n',E*10,nu,0); %Modulus - Poisson - Initial Temp
fprintf(fid,'*EXPANSION, TYPE=ISO\n');
Alpha = 0.01; % 1 percent for 1 degree
fprintf(fid,'%20.10E,%20.10E\n',Alpha,0); %Exapnsion Coeffiction - Initial Temp
fprintf(fid,'*SOLID SECTION, ELSET=incl_elements, MATERIAL=ThermalParticle\n'); 

for ijjk = 1 : size(part_nodes_idx, 1)
fprintf(fid,'*INITIAL CONDITIONS, TYPE=TEMPERATURE\n');
fprintf(fid,['CNODES_', num2str(ijjk),', 0\n']);
end


%---Load Steps---
fprintf(fid,'**\n');
fprintf(fid,'** Load Steps\n');
fprintf(fid,'**\n');

radial_strain = radial_strain / NN ;

ik2 = 0 ; % initialize

for MM = 1 : cycle

% --- STEP 1-NN ---
for ik = 1:NN
    
fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=300000\n', (MM-1)*(ik+ik2)+ik);
fprintf(fid,'*DYNAMIC, APPLICATION=QUASI-STATIC\n');
fprintf(fid,'%10.20f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 

% COMMENT IF FREE EXTERNAL BOUNDARY 
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YPOS, 1, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YPOS, 2, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YNEG, 1, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YNEG, 2, , %8.4f\n', 0); 

    % random search of inclusion:
    if fraction ~= 0
    pool = (1:1:size(x0, 1)) ;
    pqr = floor(fraction * size(x0, 1)) ;
    rst = floor(size(x0, 1)/pqr) ;
    kk = 0 ; 
    for jk = 1 : rst : size(x0, 1)
        if jk <= (size(x0, 1)-rst)
        kk = kk + 1 ;
        pool_tosel = pool((kk-1)*rst+1:(kk)*rst) ;
        selected_incl(kk) = randsample(pool_tosel,1) ;
        clear pool_tosel
        end
    end
    end
    % end of random search of inclusion...
    
    if fraction ~= 1 && fraction ~= 0
    for lmn = selected_incl 
    for ijk = ((lmn-1) * el_per_incl + 1) : ((lmn) * el_per_incl) 
    fprintf(fid,'*TEMPERATURE\n');
    fprintf(fid,['CNODES_', num2str(ijk),', %8.10f\n'], -radial_strain*(ik));
    end
    end
    elseif fraction == 1
    for ijk = 1 : size(part_nodes_idx, 1)
    fprintf(fid,'*TEMPERATURE\n');
    fprintf(fid,['CNODES_', num2str(ijk),', %8.10f\n'], -radial_strain*(ik));
    end
    end    
    
    fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=20 \n');
    fprintf(fid,'*NODE OUTPUT\n');
    fprintf(fid,'COORD, U, RF\n');
    fprintf(fid,'*ELEMENT OUTPUT\n');
    fprintf(fid,'SE, LE, S, SK, SF, NFORC\n'); % SK gives section curvature change
    fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL, NUMBER INTERVAL=20\n');

    fprintf(fid,'*END STEP\n');
    
end


% --- STEP NN+1:2NN ---
for ik2 = 1:1
    
fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=300000\n', (MM-1)*(ik+ik2) + (ik+ik2));
fprintf(fid,'*DYNAMIC, APPLICATION=QUASI-STATIC\n');
fprintf(fid,'%10.20f,%8.0f,%20.3e,%8.2f\n',[dt duration*NN 0 0]); 


% COMMENT IF FREE EXTERNAL BOUNDARY 
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YPOS, 1, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YPOS, 2, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YNEG, 1, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YNEG, 2, , %8.4f\n', 0); 

    % random search of inclusion:
    if fraction ~= 0
    pool = (1:1:size(x0, 1)) ;
    pqr = floor(fraction * size(x0, 1)) ;
    rst = floor(size(x0, 1)/pqr) ;
    kk = 0 ; 
    for jk = 1 : rst : size(x0, 1)
        if jk <= (size(x0, 1)-rst)
        kk = kk + 1 ;
        pool_tosel = pool((kk-1)*rst+1:(kk)*rst) ;
        selected_incl(kk) = randsample(pool_tosel,1) ;
        clear pool_tosel
        end
    end
    end
    % end of random search of inclusion...
    
    if fraction ~= 1 && fraction ~= 0
    for lmn = selected_incl 
    for ijk = ((lmn-1) * el_per_incl + 1) : ((lmn) * el_per_incl) 
    fprintf(fid,'*TEMPERATURE\n');
fprintf(fid,['CNODES_', num2str(ijk),', %8.10f\n'], 0);
    end
    end
    elseif fraction == 1
    for ijk = 1 : size(part_nodes_idx, 1)
    fprintf(fid,'*TEMPERATURE\n');
fprintf(fid,['CNODES_', num2str(ijk),', %8.10f\n'], 0);
    end
    end    
    
    fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=20 \n');
    fprintf(fid,'*NODE OUTPUT\n');
    fprintf(fid,'COORD, U, RF\n');
    fprintf(fid,'*ELEMENT OUTPUT\n');
    fprintf(fid,'SE, LE, S, SK, SF, NFORC\n'); % SK gives section curvature change
    fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL, NUMBER INTERVAL=20\n');

    fprintf(fid,'*END STEP\n');
    
end

end

fclose('all');

return
