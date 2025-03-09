function Aave = ave_circ(A)
%
% Average data around a circle about the center of the image
% IN FUTURE ADD CAPABILITY TO AVERAGE ABOUT AN ARBITRARY POINT
%
% INPUTS
% A     2D array of data to average
%       Requires equal grid point spacing in x and y directions
%       Assumes A has an odd number of points; otherwise center position
%       chosen will be a half gridpoint off
%
% OUTPUTS
% Aave     Averaged data, vector with same spacing as gridpoints of x
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2015-2017


% Gridpoints on which we have data
[M, N] = size(A);
[X, Y] = meshgrid(1:N,1:M);
% Square of distance from center of A
R2 = (X- (N+1)/2).^2 + (Y- (M+1)/2).^2;
% Preallocate averaged array
rmax = floor(min([M, N])/2)-1; % max radius of averaging circle
A_ave = zeros(1,rmax)*nan;
% Compute averages for each d0-spaced radius
for r=1:rmax
    % Get contour indices a distance r^2 from the origin using contourc command
    c = contourc(1:N, 1:M, R2, [0 0]+r^2);
    % First column specifies contour levels (r^2) and number of pixels
    % at that level so remove it. Also round to nearest pixels.
    c = round(c(:,2:end));
    % Get values on contour
    A_c = A(sub2ind([M,N],c(2,:),c(1,:)));
    % Average values on contour
    A_ave(r) = mean(A_c) ;
end
% Add a value for zero radius
Aave = [A((M+1)/2,(N+1)/2), A_ave];  


