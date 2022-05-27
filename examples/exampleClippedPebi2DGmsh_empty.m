%% Example of pebiGrid2DGmsh usage
% Create a simple 2D mesh using Gmsh, called from MATLAB

% Setup required variables
resGridSize = 0.1;
domain = [1 1];

% pebiGrid2DGmsh returns a single variable, G
G = clippedPebi2DGmsh(resGridSize, domain);

% Plot result
hold on;
plotGrid(G, 'faceColor', 'none');

% Save plot
f = gcf;
exportgraphics(f,'exampleClippedPebi2DGmsh_empty.png')