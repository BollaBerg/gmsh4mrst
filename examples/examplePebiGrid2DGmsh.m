%% Example of pebiGrid2DGmsh usage
% Create a simple 2D mesh using Gmsh, called from MATLAB
clear all;
% Setup required variables
resGridSize = 0.2;
domain = [
    0 0; 0.5 0.2; 1 0; 1 1; 0 1;
];
faceConstraints = {
    [0.25 0.25; 0.4 0.5; 0.7 0.7;];
    [0.8 0.1; 0.9 0.2];
    [0.2 0.9; 0.9 0.1];
    [0.1 0.1];
};

cellConstraints.a.x = 0.1;
cellConstraints.a.y = 0.1;
cellConstraints.b.x = [0.5 0.6 0.7 0.9];
cellConstraints.b.y = [0.8 0.7 0.8 0.6];

% pebiGrid2DGmsh returns a single variable, G
G = pebiGrid2DGmsh( ...
    resGridSize, ...
    domain, ...
    'faceConstraints', faceConstraints, ...
    'cellConstraints', cellConstraints, ...
    'faceConstraintFactor', 1/3, ...
    'meshAlgorithm', 'DelQuad' ...
);

% Plot result
hold on;
plotGrid(G, 'faceColor', 'none');

% Plot lines
plot(faceConstraints{1}(:, 1), faceConstraints{1}(:, 2), 'lineWidth', 2, 'color', 'magenta', 'markersize', 3);
plot(faceConstraints{2}(:, 1), faceConstraints{2}(:, 2), 'lineWidth', 2, 'color', 'magenta', 'markersize', 3);
plot(faceConstraints{3}(:, 1), faceConstraints{3}(:, 2), 'lineWidth', 2, 'color', 'magenta', 'markersize', 3);
plot(cellConstraints.a.x, cellConstraints.a.y, 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--', Marker='o');
plot(cellConstraints.b.x, cellConstraints.b.y, 'LineWidth', 1, 'Color', 'blue', 'LineStyle','--');

% Save plot
f = gcf;
exportgraphics(f,'examplePebiGrid2DGmsh.png')