clear all 
close all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDF configuration
x_min = 20; x_max = 30;
y_min = 0; y_max = 10;
e_mu_x = 10; 
e_mu_y = 10;
exp_offset = -10;
sampleSize = 1500; % total samples per distribution
HalfPeriod_L = [ x_max+(x_max - x_min)*5 y_max+(y_max - y_min)*5 ];
Mixture_N = 25; % total number of mixture components
Trimming_Factor = 5; 
Convergence_Factor = 0.000001;
% view preferences
my_xlim = 50; 
my_ylim = 50; 
my_zlim = 0.009;
my_view = [429.78244 20.99755];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
% Exponential samples
x_samples_exponential = exprnd(e_mu_x, sampleSize, 1) + exp_offset;
y_samples_exponential = exprnd(e_mu_y, sampleSize, 1) + exp_offset;
% Input concatenation & shuffle
x_input = [ x_samples_exponential ];
y_input = [ y_samples_exponential ] ;
inputArray = [x_input(:), y_input(:)];
idx = randperm(size(inputArray, 2));
inputArray = inputArray(:, idx);
% plot input in scattered form 
figure;
scatter(x_input, y_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMoFS Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pdf_index, Estimated_PDF] = emofs( ...
    inputArray, HalfPeriod_L, Mixture_N, Trimming_Factor, ...
    Convergence_Factor); 
% Generate X,Y,Z indices from pdf_index
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
surf(X, Y, Z_emofs); 
xlabel('x (in)');
ylabel('y (in)');
zlabel('Estimated PDF by EMoFS');
set(gca,'ydir','reverse');
% Show in 2D contour
figure;
contour(X, Y, Z_emofs, Mixture_N, 'k');  % 2D contour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original PDF and its illustration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential
Z_orig_exp = exppdf(X - exp_offset, e_mu_x) .* exppdf(Y - exp_offset, e_mu_y);
% Concatenation & Normalization - not necessary 
Z_orig = Z_orig_exp;
% Plot
figure;
surf(X, Y, Z_orig); 
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF');
set(gca,'ydir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw EMoFS and Original Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
axEMOFS = ax;
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
xlim([ -my_xlim my_xlim ]);
ylim([ -my_ylim my_ylim ]);
zlim([0 my_zlim]);
view( my_view );
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF vs EMoFS');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IGMM Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = statset('MaxIter',1000, 'Display','final');
gmdist = fitgmdist(inputArray, Mixture_N,'Options',options);
grid_points = [x_axis, y_axis];
f = pdf(gmdist, grid_points);
Z_igmm = reshape(f, MaxPeriod, MaxPeriod);
% Plot
figure;
h1 = surf(X, Y, Z_igmm);
xlabel('x (in)');
ylabel('y (in)');
zlabel('Estimated PDF by IGMM');
set(gca,'ydir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw EMoFS and Original Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
axEMOFS = ax;
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
xlim([ -my_xlim my_xlim ]);
ylim([ -my_ylim my_ylim ]);
zlim([0 my_zlim]);
view( my_view );
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF vs EMoFS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw IGMM and Original Together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the first surface using colors from the winter colormap
fig = figure(); 
tcl = tiledlayout(fig,1,1); 
ax = nexttile(tcl); 
axIGMM = ax;
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
xlim([ -my_xlim my_xlim ]);
ylim([ -my_ylim my_ylim ]);
zlim([0 my_zlim]);
view( my_view );
xlabel('x (in)');
ylabel('y (in)');
zlabel('Original PDF vs IGMM');

%linkprop([axEMOFS, axIGMM], ...
%    {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMS calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_diff = Z_orig - Z_emofs;
error_rms = sqrt(sum(error_diff.^2, "all")); 
fprintf("RMS EMoFS:  %f\n", error_rms);
fprintf("Total Sum EMoFS: %f\n", sum(Z_emofs, "all"));
error_diff = Z_orig - Z_igmm;
error_rms = sqrt(sum(error_diff.^2, "all")); 
fprintf("RMS IGMM:  %f\n", error_rms);
fprintf("Total Sum IGMM: %f\n", sum(Z_igmm, "all"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormapToTruecolor description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tc = colormapToTruecolor(map,ZData)
    % map is a n-by-3 colormap matrix
    % ZData is a k-by-w matrix of ZData (or CData, I suppose)
    % tc is a kxwx3 Truecolor array based on map and ZData values.
    tcIdx = round(rescale(ZData,1,height(map)));
    tc = reshape(map(tcIdx,:),[size(ZData),3]);
end