function [f_FT, U_FT] = nufft2d_universal_ALTnorm(U, x, y, Nx, Ny)
%nufft2d_universal
% 
% [f_FT, U_FT] = nufft2d_universal(U, x, y, Nx, Ny)
% 
% Compute a 2D fourier transform on a nonuniform set of data points in rectangular domain. 
% 
% IMPORTANT: This function works for any rectangular domain (hence named 'universal')
% 
% This code normalizes the Fourier transformed data by the zero frequncy
% data, which is proportional to the mean of u
%
% INPUTS 
% U    data of interest (column vector)
% x    x data points (column vector of same length as U)
% y    y data points (column vector of same length as U)
% Nx   Number of data points along x on which to compute fourier transform. Must be 
%      odd. Also should match your spatial resolution. 
% Ny   Number of data points along y on which to compute fourier transform. Must be 
%      odd. Also should match your spatial resolution. 
% 
% OUTPUTS
% f   vector of spatial frequencies. Units: inverse of units of x and y.
% U   Fourier transformed data. Units: None, because data was normalized as
%     described above. Length of this vector is determined by minimum of Nx
%     and Ny. 
% 
% Requires the function ave_circ.m
% 
% Mainak Sarkar, Aug 2, 2021 

% --- PREPARE INPUTS ----

% Scale x and y data to go from 0 to 1 -- seems to be required based on
% reading comments in the m file nufftn
xscale = x/max(x);
yscale = y/max(y);

% --- FFT WITH CUSTOM QUERY POINTS ---

% There are N sample points in the X and Y directions. Hence sampling
% periods are
Tx = (max(x)-min(x))/Nx ; % units: um (x and y are not scaled to go 0 to 1)
Ty = (max(y)-min(y))/Ny ; 

% Query points in both x and y 
Qx = 0 : 1 : (Nx-1) ;
Qy = 0 : 1 : (Ny-1) ;
[Qg1, Qg2] = meshgrid(Qx, Qy) ;

ctr = 0 ;
for cm = 1 : size(Qg1, 1)
    for cn = 1 : size(Qg1, 2)
        ctr = ctr + 1 ;
        query_pts(ctr, 1:2) = [ Qg1(cm, cn), Qg2(cm, cn) ] ;
    end
end

% Nonuniform fft
U_FT0 = nufftn(U, [xscale, yscale], query_pts) ;

% --- RESHAPE FFT DATA SET ---
% U_FT0 is a vector corresponding to a 2D data set. Reshape it into a 2D
% array

for jj = 1 : size(query_pts, 1)
    mm = query_pts(jj, 1) + 1 ;
    nn = query_pts(jj, 2) + 1 ;
    U_FT2D(nn, mm) = U_FT0(jj) ;
end

U_FT2D = abs(fftshift(U_FT2D)) ;

%{
% Scale so max value is 1
U_FT2D = U_FT2D/(mean(U)*numel(U));
%}
%{
% Alternative normalization to make energy analogous to probability: 
U_FT2D = U_FT2D./sum(U_FT2D(:)) ; 
%}
% Alternative: Do not normalize: 
U_FT2D = U_FT2D ; 

% use contour plot for debugging (if necessary)
% contour(U_FT2D)

% Average data around circle to get a vector
U_FT = ave_circ(U_FT2D) ;

% --- Vector defining sampling frequency ---
f_FT = (0:(min([Nx Ny])-1)/2-1)/(Ty*min([Nx Ny])) ;

return

