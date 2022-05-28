%% Example of pebiGrid2DGmsh usage
% Create a simple 2D mesh using Gmsh, called from MATLAB
clear all;

% Setup required variables
resGridSize = 0.1;
domain = [
    0 0; 0.5 0.2; 1 0; 1 1; 0 1;
];
faceConstraints = {
    [0.25 0.25; 0.4 0.5; 0.7 0.7;], ...
    [0.8 0.1; 0.9 0.2], ...
};
cellConstraints = {
    [0.1 0.1], ...
    [0.5 0.8; 0.6 0.7; 0.7 0.8; 0.9 0.6];
};

% pebiGrid2DGmsh returns a single variable, G
G = pebiGrid2DGmshBase( ...
    resGridSize, ...
    domain, ...
    'faceConstraints', faceConstraints, ...
    'cellConstraints', cellConstraints ...
);

% Plot result
hold on;
plotGrid(G, 'faceColor', 'none');

% Plot lines
plot(faceConstraints{1}(:, 1), faceConstraints{1}(:, 2), 'lineWidth', 2, 'color', 'magenta', 'markersize', 3, Marker='+', LineStyle='none');
plot(faceConstraints{2}(:, 1), faceConstraints{2}(:, 2), 'lineWidth', 2, 'color', 'magenta', 'markersize', 3, Marker='*', LineStyle='none');
plot(cellConstraints{1}(:, 1), cellConstraints{1}(:, 2), 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--', Marker='o', LineStyle='none');
plot(cellConstraints{2}(:, 1), cellConstraints{2}(:, 2), 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--', Marker='.', LineStyle='none');

% Save plot
f = gcf;
exportgraphics(f,'examplePebiGrid2DGmshBase.png')