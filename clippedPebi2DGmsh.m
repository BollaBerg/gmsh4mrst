function G = clippedPebi2DGmsh(resGridSize, size, varargin)
% Construct a 2D clipped PEBI grid, using Gmsh
% 
% SYNOPSIS:
%   G = pebiGrid2DGmsh(resGridSize, pdims, faceConstraints)
%
% ARGUMENTS
%   resGridSize     - Size of the reservoir grid cells, in units of meters.
%
%   size            - Vector, length 2, [xmax, ymax], of physical size in
%                   units of meters of the computational domain. The
%                   domain always starts at [0, 0].
%
% OPTIONAL PARAMETERS
%   faceConstraints - A struct of vectors. Each vector, size nf x 2, is the
%                   coordinates of a surface-trace. The surface is
%                   assumed to be linear between the coordinates. The
%                   function will place sites such that the surface is
%                   traced by faces of the grid. Defaults to empty
%                   struct.
%
%   faceConstraintFactor - Float. The size of the cells close to the face
%                   constraints, as compared to supplied cell_dimensions.
%                   Cells within min_threshold_distance will have size
%                   face_constraint_factor * cell_dimensions. Equivalent
%                   to FCFactor in MRST/UPR/pebiGrid2D. Defaults to 1/3.
%
%   minThresholdDistance - Float. Distance from face constraints where cell
%                   dimensions will start increasing. Defaults to 0.05.
%
%   maxThresholdDistance - Float. Distance from face constraints where cell
%                   dimensions will be back to their default (max) value,
%                   i.e. the supplied argument cell_dimensions. Defaults
%                   to 0.2.
%
%   faceIntersectionFactor - Float. The size of the cells close to
%                   intersections between face constraints, as compared to supplied
%                   cell_dimensions. Cells in min_intersection_distance
%                   from an intersection will have size
%                   face_intersection_factor * cell_dimensions.
%                   If missing, no extra cell shaping will occur around
%                   intersections. The factor is also used in "breaks" of
%                   lines, i.e. if there is a sharp "turn" in a line
%                   segment. Defaults to missing.
%
%   minIntersectionDistance - Float. Distance from intersections where cell
%                   dimensions will start increasing. If missing, will use
%                   min_threshold_distance. Defaults to missing.
%
%   maxIntersectionDistance - Float. Distance from intersections where
%                   cell dimensions will be back to their default (max)
%                   value, i.e. the supplied argument cell_dimensions. If
%                   missing, will use min_threshold_distance. Defaults to
%                   missing.
%
%   fractureMeshSampling - Int. The number of points along the face
%                   constraints should be sampled to calculate the
%                   threshold distances. Defaults to 100.
%
%   cellConstraints - A struct of vectors. Each vector, size nf x 2, is a
%                   line, which the method will attempt to place in the
%                   center of the returned grid. The lines are assumed to
%                   be linear between the coordinates. If a constraint is
%                   only one coordinate, the line is treated as a point
%                   constraint. Defaults to empty struct.
%
%   cellConstraintFactor - Float. The size used for cells around cell
%                   constraints, as compared to supplied cell_dimensions.
%                   Equivalent to CCFactor in MRST/UPR/pebiGrid2D.
%                   Defaults to 1/4.
%
%   cellConstraintLineFactor - Float. The size used for cells around cell
%                   constraint lines, as compared to the supplied
%                   cell_dimensions. Overrides cell_constraint_factor for
%                   lines. If missing, cell_constraint_factor will be
%                   used for lines. Defaults to missing.
%
%   cellConstraintPointFactor - Float. The size used for cells around cell
%                   constraint points, as compared to the supplied 
%                   cell_dimensions. Overrides cell_constraint_factor for
%                   points. If missing, cell_constraint_factor will be used
%                   for points. Defaults to missing.
%
%   meshAlgoritm - String | Int. What meshing algorithm should be used.
%                   Can either be the Gmsh-given ID of the algorithm or the
%                   name of the algorithm. Legal values:
%                       "MeshAdapt" = 1
%                       "Automatic" = 2
%                       "Delaunay" = 5
%                       "Frontal" = 6
%                       "BAMG" = 7
%                       "DelQuad" = 8
%                   Defaults to "Delaunay".
%
%   recombinationAlgorithm - String | Int. What recombination algorithm
%                   should be used, and whether recombination should be
%                   done. Recombination makes Gmsh attempt to create a
%                   quadrangle mesh, rather than a triangle mesh. Can be
%                   either the Gmsh-given ID of the algorithm or the name
%                   of the algorithm. Legal values:
%                       "Simple" = 0
%                       "Blossom" = 1
%                       "SimpleFull" = 2
%                       "BlossomFull" = 3
%                   If missing, no recombination will be done. Defaults to
%                   missing.
%
%                   NOTE: Blossom is Gmsh default, but struggles with
%                       constraints. It may therefore be beneficial to use
%                       SimpleQuad, if there are constraints.
%
%                   NOTE 2: Recombination may lead to weird results for
%                       constraints, as the recombination is done after
%                       constraints have been applied. This is especially
%                       true for SimpleFull and BlossomFull, who 
%                       automatically perform a coarser mesh, followed by
%                       recombination, smoothing and subdivision.
%


% clippedPebi2D assumes P is an array of Voronoi sites, while 
% pyG.nodes.coords is an array of Delaunay sites. However, due to the
% duality of the two, we can use pyG.nodes.coords as Voronoi sites,
% as long as we swap faceConstraints and cellConstraints. This way, the
% Delaunay sites (in pyG.nodes.coords) are located along the
% faceConstraints, which gives us the right result when using them as
% Voronoi sites in clippedPebi2D.
for i = 1:length(varargin)
    if strcmp(varargin(i), 'faceConstraints')
        varargin{i} = 'cellConstraints';
    elseif strcmp(varargin(i), 'cellConstraints')
        varargin{i} = 'faceConstraints';
    elseif strcmp(varargin(i), 'cellConstraintFactor')
        varargin{i} = 'faceConstraintFactor';
    elseif strcmp(varargin(i), 'faceConstraintFactor')
        varargin{i} = 'cellConstraintFactor';
    end
end

pyG = pebiGrid2DGmsh(resGridSize, size, varargin{:});

bnd = [
    0 0;
    size(1) 0;
    size(1) size(2);
    0 size(2);
];

G = clippedPebi2D(pyG.nodes.coords, bnd);

end