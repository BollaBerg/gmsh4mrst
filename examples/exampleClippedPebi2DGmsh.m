%% Example of pebiGrid2DGmsh usage
% Create a simple 2D mesh using Gmsh, called from MATLAB
clear all;
% Setup required variables
resGridSize = 0.2;
domain = [1 1];
faceConstraints.x = [0.1 0.5 0.9];
faceConstraints.y = [0.8 0.6 0.2];

cellConstraints.a.x = 0.1;
cellConstraints.a.y = 0.1;
cellConstraints.b.x = [0.5 0.6 0.7 0.9];
cellConstraints.b.y = [0.8 0.7 0.8 0.6];

% pebiGrid2DGmsh returns a single variable, G
G = clippedPebi2DGmsh( ...
    resGridSize, ...
    domain, ...
    'faceConstraints', faceConstraints, ...
    'cellConstraints', cellConstraints ...
);

% Plot result
hold on;
plotGrid(G, 'faceColor', 'none');

% Plot lines
plot(faceConstraints.x, faceConstraints.y, 'lineWidth', 2, 'color', 'magenta', 'markersize', 3);
plot(cellConstraints.a.x, cellConstraints.a.y, 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--', Marker='o');
plot(cellConstraints.b.x, cellConstraints.b.y, 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--');

% Save plot
f = gcf;
exportgraphics(f,'exampleClippedPebi2DGmsh.png')