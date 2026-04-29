clear all 
close all 

% view settings
my_xlim = 60; 
my_ylim = 60; 
my_zlim = 0.011;
my_view = [324 15];

% Define 2D uniform range
x_min = 10; x_max = 20;
y_min = 10; y_max = 20; 

% Define parameters
sampleSize = 1500;
HalfPeriod_L = [ x_max+(x_max - x_min)*3 y_max+(y_max - y_min)*3 ];
Mixture_N = 10;
Trimming_Factor = 20; 
Convergence_Factor = 0.000001;

% Gaussian parameters
mu_x = (x_min + x_max) / 2;
mu_y = (y_min + y_max) / 2;
sigma_x = (x_max - x_min) / 2;
sigma_y = (y_max - y_min) / 2;
% Generate Gaussian samples
x_samples = mu_x + sigma_x * randn(sampleSize, 1);
y_samples = mu_y + sigma_y * randn(sampleSize, 1);

% Finalize input
inputArray = [x_samples y_samples];

% plot output
figure;
scatter(x_samples, y_samples);

% estimate the pdf with EMoFS
[pdf_index, Estimated_PDF] = emofs( ...
    inputArray, HalfPeriod_L, Mixture_N, Trimming_Factor, ...
    Convergence_Factor); 

% Prepare X,Y,Z indices 
MaxPeriod = max(HalfPeriod_L) * 2 + 1;
x_axis = pdf_index(:,1);
y_axis = pdf_index(:,2);
z_axis = Estimated_PDF;
% Reshape into grid form
X = reshape(x_axis, MaxPeriod, MaxPeriod);
Y = reshape(y_axis, MaxPeriod, MaxPeriod);
Z_emofs = reshape(z_axis, MaxPeriod, MaxPeriod);
% Plot the estimated pdf
figure;
%surf(X, Y, Z, ...
%    'FaceColor', 'none', ...   % remove surface fill
%    'EdgeColor', 'interp');    % color edges based on Z
surf(X, Y, Z_emofs); 
xlabel('x (ph)');
ylabel('y (ph)');
zlabel('Estimated PDF by EMoFS');
set(gca,'ydir','reverse');

% calculate original PDF and draw it
%x_orig = linspace(mu_x - 3*sigma_x, mu_x + 3*sigma_x, 100);
%y_orig = linspace(mu_y - 3*sigma_y, mu_y + 3*sigma_y, 100);
%[X_orig, Y_orig] = meshgrid(x_orig, y_orig);
Z_orig = (1 / (2*pi*sigma_x*sigma_y)) * ...
      exp(-((X - mu_x).^2 / (2*sigma_x^2)+...
      (Y - mu_y).^2 / (2*sigma_y^2)));
figure;
surf(X, Y, Z_orig); 
xlabel('x (ph)');
ylabel('y (ph)');
zlabel('Original PDF');
set(gca,'ydir','reverse');

% 2D contour
figure;
contour(X, Y, Z_emofs, Mixture_N, 'k');  % 2D contour

% estimate the pdf with IGMM
options = statset('MaxIter',1000, 'Display','final');
gmdist = fitgmdist(inputArray, Mixture_N,'Options',options); % Mixture_N components
grid_points = [x_axis, y_axis];
f = pdf(gmdist, grid_points);
Z_igmm = reshape(f, MaxPeriod, MaxPeriod);

figure;
h1 = surf(X, Y, Z_igmm);
%set(h1, 'FaceAlpha', 0.6);
xlabel('x (ph)');
ylabel('y (ph)');
zlabel('Estimated PDF by IGMM');
set(gca,'ydir','reverse');
%legend('GMM', 'EMoFS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw EMoFS and IGMM Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
surf1Colormap = winter(); 
tc = colormapToTruecolor(surf1Colormap,Z_igmm);  
hsurf1 = surf(ax,X,Y,Z_igmm,tc,'FaceColor','interp','EdgeAlpha',0.3);
set(hsurf1, 'FaceAlpha', 0.2);
% create second surface using colors from the autumn colormap
hold on
surf2Colormap = spring(); 
tc = colormapToTruecolor(surf2Colormap,Z_emofs);
hsurf2 = surf(ax,X,Y,Z_emofs,tc,'FaceColor','interp','EdgeAlpha',0.3);
axHidden = axes(tcl,'visible','off','HandleVisibility','off');
colormap(ax,surf1Colormap)
colormap(axHidden,surf2Colormap)
cb1 = colorbar(ax);
cb1.Layout.Tile = 'east';
cb1.Label.String = 'Estimation by IGMM';
cb2 = colorbar(axHidden);
cb2.Layout.Tile = 'east';
cb2.Label.String = 'Estimation by EMoFS';
xlim([ -(x_min*2) (x_max * 2) ]);
ylim([ -(y_min*2) (y_max * 2) ]);
zlim([0 max(Z_igmm(:))]);
xlabel('x (in)');
ylabel('y (in)');
zlabel('IGMM vs EMoFS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw EMoFS and Original Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
surf1Colormap = winter(); 
tc = colormapToTruecolor(surf1Colormap,Z_orig);  
hsurf1 = surf(ax,X,Y,Z_orig,tc,'FaceColor','interp','EdgeAlpha',0.3);
set(hsurf1, 'FaceAlpha', 0.2);
% create second surface using colors from the autumn colormap
hold on
surf2Colormap = spring(); 
tc = colormapToTruecolor(surf2Colormap,Z_emofs);
hsurf2 = surf(ax,X,Y,Z_emofs,tc,'FaceColor','interp','EdgeAlpha',0.3);
axHidden = axes(tcl,'visible','off','HandleVisibility','off');
colormap(ax,surf1Colormap)
colormap(axHidden,surf2Colormap)
cb1 = colorbar(ax);
cb1.Layout.Tile = 'east';
cb1.Label.String = 'Original PDF';
cb2 = colorbar(axHidden);
cb2.Layout.Tile = 'east';
cb2.Label.String = 'Estimation by EMoFS';
xlim([ -(x_min*2) (x_max * 2) ]);
ylim([ -(y_min*2) (y_max * 2) ]);
zlim([0 max(Z_igmm(:))]);
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF vs EMoFS');
view( my_view );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw IGMM and Original Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
surf1Colormap = winter(); 
tc = colormapToTruecolor(surf1Colormap,Z_orig);  
hsurf1 = surf(ax,X,Y,Z_orig,tc,'FaceColor','interp','EdgeAlpha',0.3);
set(hsurf1, 'FaceAlpha', 0.2);
% create second surface using colors from the autumn colormap
hold on
surf2Colormap = spring(); 
tc = colormapToTruecolor(surf2Colormap,Z_igmm);
hsurf2 = surf(ax,X,Y,Z_igmm,tc,'FaceColor','interp','EdgeAlpha',0.3);
axHidden = axes(tcl,'visible','off','HandleVisibility','off');
colormap(ax,surf1Colormap)
colormap(axHidden,surf2Colormap)
cb1 = colorbar(ax);
cb1.Layout.Tile = 'east';
cb1.Label.String = 'Original PDF';
cb2 = colorbar(axHidden);
cb2.Layout.Tile = 'east';
cb2.Label.String = 'Estimation by IGMM';
xlim([ -(x_min*2) (x_max * 2) ]);
ylim([ -(y_min*2) (y_max * 2) ]);
zlim([0 max(Z_igmm(:))]);
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF vs IGMM');
view( my_view );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMS calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_diff = Z_orig - Z_emofs;
error_rms = sqrt(sum(error_diff.^2, "all")); 
fprintf("RSS EMoFS:  %f\n", error_rms);
fprintf("Total Sum EMoFS: %f\n", sum(Z_emofs, "all"));
error_diff = Z_orig - Z_igmm;
error_rms = sqrt(sum(error_diff.^2, "all")); 
fprintf("RSS IGMM:  %f\n", error_rms);
fprintf("Total Sum IGMM: %f\n", sum(Z_igmm, "all"));


function tc = colormapToTruecolor(map,ZData)
% map is a n-by-3 colormap matrix
% ZData is a k-by-w matrix of ZData (or CData, I suppose)
% tc is a kxwx3 Truecolor array based on map and ZData values.
tcIdx = round(rescale(ZData,1,height(map)));
tc = reshape(map(tcIdx,:),[size(ZData),3]);
end