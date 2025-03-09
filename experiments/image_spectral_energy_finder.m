% This script will measure the total spectral energy of the inputted matrix
% image. 
% Mainak Sarkar, 2024

clear all
clc
close all

% Add your confocal microscopic image of matrix:
I = imread('C:\...\sample_matrix_image.tif') ;
size(I)


subset_size = 32-1 ; % px; keep odd.

count = 0 ;
countx = 0 ;
for i = 1 : subset_size : size(I, 2)-(subset_size-1)
    countx = countx + 1 ;
    county = 0 ;
    for j = 1 : subset_size : size(I, 1)-(subset_size-1)
        county = county + 1 ;
        count = count + 1 ;
        Is{count} = I(j:j+subset_size-1, i:i+subset_size-1) ;
        Is_double{count} = double(Is{count});
        Is_double{count} = rescale(Is_double{count}) ;
        id_subimg(count,:) = [count countx county] ;
    end
end

% check point
% imshow(Is{1000})

 
for lmn = 1 : size(Is_double, 2)
% FFT:
[xG,yG]=meshgrid(1:size(Is_double{lmn},2), 1:size(Is_double{lmn},1));


% FID
xG = xG - min(xG(:)) ;
yG = yG - min(yG(:)) ;
u_FID = Is_double{lmn} ;

% return
% xG = xG(2:end-1, 2:end-1) ;
% yG = yG(2:end-1, 2:end-1) ;
% u_FID = u_FID(2:end-1, 2:end-1) ;

% Keep non nan data
I=~isnan(u_FID);
XX=xG(I); YY=yG(I); UU=u_FID(I);

% choose resolution (must be odd)
Nx = size(xG, 2) ;
Ny = size(xG, 1) ;

% return
[f_FT_fidCH{lmn}, U_FT_fidCH{lmn}] = nufft2d_universal_ALTnorm(UU, XX, YY, Nx, Ny);

% measure power and energy:
n = length(f_FT_fidCH{lmn});          % number of samples
power{lmn} = abs(U_FT_fidCH{lmn}).^2/n;    % power of the DFT
TotalEnergy(lmn) = sum(power{lmn}(:)) ;

end



sum(TotalEnergy(:))
std(TotalEnergy)
mean(TotalEnergy)
mean(TotalEnergy)/std(TotalEnergy)


% plot:
a2 = histogram(TotalEnergy) ; 
b2 = a2.Values ;
c2 = a2.BinEdges ;
close
hf = make_fig([1 1 0.5 0.5]) ;
plot(c2(1:end-1), b2, '-','linewidth',1.5,'Color',[0 0 0]) ;
hold on
set(gcf,'units','normalized','outerposition',[0 0 1 1])
pbaspect([1 1 1])
ax = gca;
set(gca, 'box', 'off')
ylabel('#')
xlabel('spectral energy, subimage')
% ylim([0 750])
% xlim([0 4000])
% set(gca, 'YScale', 'log', 'xscale', 'log')
% legend('inclusions recovered','inclusions contracted')
% legend('box', 'off', 'Location', 'northwest')




for ijk = 1 : size(id_subimg,1)
    powermesh(id_subimg(ijk,3), id_subimg(ijk,2)) = TotalEnergy(id_subimg(ijk,1)) ;
end


[xGi, yGi] = meshgrid(1:size(powermesh,2), 1:size(powermesh,1)) ;

hf = make_fig([1 1 0.5 0.5]) ;
surf(xGi, yGi, powermesh)
view(2)
daspect([1 1 1]) ;
shading interp
colormap(jet(2048))
% colormap(redblue)
% colormap(parula)
h = colorbar ;
ax = gca ;
jj = ax.CLim ;
cmin = jj(1) ;
cmax = jj(2) ;
lt2 = max(abs(cmin), abs(cmax)) ;
set(gca, 'fontsize',12,'fontweight','bold')
ylabel(h, '|Energy|', 'fontsize',10,'fontweight','bold')
box off
grid off
axis off
caxis([0 lt2])



return
































