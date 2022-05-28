function G = pebiGrid2DGmshBase(resGridSize, size, varargin)
% Construct a 2D PEBI grid, using Gmsh.
%
% This method creates a Gmsh mesh using pebi_base_2D, then converts the
% mesh into a MRST PEBI grid.
% 
% SYNOPSIS:
%   G = pebiGrid2DGmshBase(resGridSize, varargin)
%   G = pebiGrid2DGmshBase(..., 'Name1', Value1, ...)
%
% ARGUMENTS
%   resGridSize     - Size of the reservoir grid cells, in units of meters.
%
%   shape           - Vector, length 2, [xmax, ymax], of physical size in
%                   units of meters of the computational domain OR
%                   - k x 2 array of coordinates. Each coordinate
%                   corresponds to a vertex in the polygon boundary. The
%                   coordinates must be ordered clockwise or counter
%                   clockwise.
%
% OPTIONAL PARAMETERS
%   faceConstraints - A cell array of float arrays. Each array, size nf x 2
%                   is the coordinates of a surface-trace. The surface is
%                   assumed to be linear between the coordinates. The
%                   function will place sites such that the surface is
%                   traced by faces of the grid. Defaults to empty cell
%                   array.
%
%   faceConstraintFactor - Float. The size of the background cells close to
%                   the face constraints, as compared to supplied
%                   resGridSize. Cells within minThresholdDistance of a
%                   face constraint will have size faceConstraintFactor *
%                   resGridSize. Only applies to the background grid.
%                   Defaults to 1/3.
%
%   faceConstraintParallelFactor - Float. The size used along cells around
%                   face constraint lines, as compared to the supplied
%                   resGridSize. Overrides faceConstraintFactor along the
%                   cells. If missing, faceConstraintFactor will be
%                   used. Defaults to missing.
%
%   faceConstraintPerpendicularFactor - Float. The size used across cells
%                   around face constraint lines, as compared to the
%                   supplied resGridSize. Overrides faceConstraintFactor
%                   across the cells. If missing, faceConstraintFactor will
%                   be used. Defaults to missing.
%
%   faceConstraintPerpendicularCells - Int. The number of transfinite cells
%                   to use across the face constraints. Should be odd if
%                   the Delaunay grid is used directly, and even if the
%                   grid is converted to a PEBI grid before use. More cells
%                   may make a "buffer zone" around the constraints.
%                   Defaults to 3.
%
%   faceConstraintPointFactor - Float. The size used for cells around face
%                   constraint points, as compared to the supplied 
%                   resGridSize. Overrides faceConstraintFactor for
%                   points. If missing, faceConstraintFactor will be used
%                   for points. Defaults to missing.
%
%   faceConstraintRefinementFactor - Float. The cell size factor that
%                   should be used for cells around the face constraint
%                   lines. If missing, resGridSize will be used and no
%                   refinement will be done. Default to missing.
%
%   minFCThresholdDistance - Float. Distance from face constraints where
%                   cell dimensions will start scaling size. Defaults to
%                   0.05.
%
%   maxFCThresholdDistance - Float. Distance from face constraints where
%                   cell dimensions will be back to their default (max)
%                   value, i.e. the supplied argument resGridSize. Defaults
%                   to 0.2.
%
%   faceMeshSampling - Int. The number of points along the face constraints
%                   should be sampled to calculate the threshold distances.
%                   Defaults to 100.
%
%   cellConstraints - A cell array of float arrays. Each array, size nf x 2
%                   is a line, which the method will attempt to place in
%                   the center of cells in the returned grid. The lines are
%                   assumed to be linear between the coordinates. If a
%                   constraint is only one coordinate, the line is treated
%                   as a point constraint. Defaults to empty cell array.
%
%   cellConstraintFactor - Float. The size used for cells around cell
%                   constraints, as compared to supplied resGridSize.
%                   Cells within minThresholdDistance of a
%                   cell constraint will have size cellConstraintFactor *
%                   resGridSize. Only applies to the background grid.
%                   Defaults to 1/4.
%
%   cellConstraintParallelFactor - Float. The size used along cells around
%                   cell constraint lines, as compared to the supplied
%                   resGridSize. Overrides cellConstraintFactor along the
%                   cells. If missing, cellConstraintFactor will be
%                   used. Defaults to missing.
%
%   cellConstraintPerpendicularFactor - Float. The size used across cells
%                   around cell constraint lines, as compared to the
%                   supplied resGridSize. Overrides cellConstraintFactor
%                   across the cells. If missing, cellConstraintFactor will
%                   be used. Defaults to missing.
%
%   cellConstraintPerpendicularCells - Int. The number of transfinite cells
%                   to use across the cell constraints. Should be even if
%                   the Delaunay grid is used directly, and odd if the grid
%                   is converted to a PEBI grid before use. More cells may
%                   make a "buffer zone" around the constraints. Defaults
%                   to 3.
%
%   cellConstraintPointFactor - Float. The size used for cells around cell
%                   constraint points, as compared to the supplied 
%                   resGridSize. Overrides cellConstraintFactor for
%                   points. If missing, cellConstraintFactor will be used
%                   for points. Defaults to missing.
%
%   cellConstraintRefinementFactor - Float. The cell size factor that
%                   should be used for cells around the cell constraint
%                   lines. If missing, resGridSize will be used and no
%                   refinement will be done. Default to missing.
%
%   minCCThresholdDistance - Float. Distance from cell constraints where
%                   cell dimensions will start scaling size. Defaults to
%                   0.05.
%
%   maxCCThresholdDistance - Float. Distance from cell constraints where
%                   cell dimensions will be back to their default (max)
%                   value, i.e. the supplied argument resGridSize. Defaults
%                   to 0.2.
%
%   CCMeshSampling - Int. The number of points along the cell constraints
%                   should be sampled to calculate the threshold distances.
%                   Defaults to 100.
%
%   convertToPEBI - Boolean. Whether the mesh should be converted to a PEBI
%                   grid before being returned. If False, the resulting
%                   grid is a Delaunay triangulation. Defaults to true.
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
%   recombinationAlgorithm - String | Int | missing. What recombination
%                   algorithm should be used, and whether recombination
%                   should be done. Recombination makes Gmsh attempt to
%                   create a quadrangle mesh, rather than a triangle mesh.
%                   Only applied to the background grid.
%                   Can be either the Gmsh-given ID of the algorithm or the
%                   name of the algorithm. Legal values:
%                       "Simple" = 0
%                       "Blossom" = 1
%                       "SimpleFull" = 2
%                       "BlossomFull" = 3
%                       missing
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

defaultFaceConstraints = {};
defaultFaceConstraintFactor = 1/3;
defaultFaceConstraintParallelFactor = string(missing);   % => Python None
defaultFaceConstraintPerpendicularFactor = string(missing);
defaultFaceConstraintPerpendicularCells = 2;
defaultFaceConstraintPointFactor = string(missing);
defaultFaceConstraintRefinementFactor = string(missing);
defaultMinFCThresholdDistance = 0.05;
defaultMaxFCThresholdDistance = 0.2;
defaultFaceMeshSampling = 100;
defaultCellConstraints = {};
defaultCellConstraintFactor = 1/4;
defaultCellConstraintParallelFactor = string(missing);
defaultCellConstraintPerpendicularFactor = string(missing);
defaultCellConstraintPerpendicularCells = 3;
defaultCellConstraintPointFactor = string(missing);
defaultCellConstraintRefinementFactor = string(missing);
defaultMinCCThresholdDistance = 0.05;
defaultMaxCCThresholdDistance = 0.2;
defaultCCMeshSampling = 100;
defaultConvertToPEBI = true;
defaultMeshAlgorithm = "Delaunay";
defaultRecombinationAlgorithm = string(missing);

    function valid = validMeshAlgorithm(x)
        legalStrings = {'MeshAdapt', 'Automatic', 'Delaunay', 'Frontal', 'BAMG', 'DelQuad'};
        if (isinteger(x) && ismember(x, [1 2 5 6 7 8]))
            valid = true;
        else
            switch validatestring(x, legalStrings)
                case legalStrings
                    valid = true;
                otherwise
                    valid = false;
            end
        end
    end

    function valid = validRecombinationAlgorithm(x)
        legalStrings = {'Simple', 'Blossom', 'SimpleFull', 'BlossomFull'};
        if ismissing(x)
            valid = true;
        elseif (isinteger(x) && ismember(x, [0 1 2 3]))
            valid = true;
        else
            switch validatestring(x, legalStrings)
                case legalStrings
                    valid = true;
                otherwise
                    valid = false;
            end
        end
    end

    validFloat = @(x) isfloat(x) && (x > 0);
    validSize = @(x) isnumeric(x) && length(x) >= 2;
    function valid = validOptionalFloat(x)
        if ismissing(x)
            valid = true;
        elseif (isfloat(x) && (x > 0))
            valid = true;
        else
            valid = false;
        end
    end
    validInt = @(x) isinteger(x) && (x > 0);

p = inputParser;
addRequired(p, 'resGridSize', validFloat);
addRequired(p, 'shape', validSize);
addParameter(p, 'faceConstraints', defaultFaceConstraints);
addParameter(p, 'faceConstraintFactor', defaultFaceConstraintFactor, validFloat);
addParameter(p, 'faceConstraintParallelFactor', defaultFaceConstraintParallelFactor, @validOptionalFloat);
addParameter(p, 'faceConstraintPerpendicularFactor', defaultFaceConstraintPerpendicularFactor, @validOptionalFloat);
addParameter(p, 'faceConstraintPerpendicularCells', defaultFaceConstraintPerpendicularCells, validInt);
addParameter(p, 'faceConstraintPointFactor', defaultFaceConstraintPointFactor, @validOptionalFloat);
addParameter(p, 'faceConstraintRefinementFactor', defaultFaceConstraintRefinementFactor, @validOptionalFloat);
addParameter(p, 'minFCThresholdDistance', defaultMinFCThresholdDistance, validFloat);
addParameter(p, 'maxFCThresholdDistance', defaultMaxFCThresholdDistance, validFloat);
addParameter(p, 'faceMeshSampling', defaultFaceMeshSampling, validInt);
addParameter(p, 'cellConstraints', defaultCellConstraints);
addParameter(p, 'cellConstraintFactor', defaultCellConstraintFactor, validFloat);
addParameter(p, 'cellConstraintParallelFactor', defaultCellConstraintParallelFactor, @validOptionalFloat);
addParameter(p, 'cellConstraintPerpendicularFactor', defaultCellConstraintPerpendicularFactor, @validOptionalFloat);
addParameter(p, 'cellConstraintPerpendicularCells', defaultCellConstraintPerpendicularCells, validInt);
addParameter(p, 'cellConstraintPointFactor', defaultCellConstraintPointFactor, @validOptionalFloat);
addParameter(p, 'cellConstraintRefinementFactor', defaultCellConstraintRefinementFactor, @validOptionalFloat);
addParameter(p, 'minCCThresholdDistance', defaultMinCCThresholdDistance, validFloat);
addParameter(p, 'maxCCThresholdDistance', defaultMaxCCThresholdDistance, validFloat);
addParameter(p, 'CCMeshSampling', defaultCCMeshSampling, validInt);
addParameter(p, 'convertToPEBI', defaultConvertToPEBI, @islogical);
addParameter(p, 'meshAlgorithm', defaultMeshAlgorithm, @validMeshAlgorithm);
addParameter(p, 'recombinationAlgorithm', defaultRecombinationAlgorithm, @validRecombinationAlgorithm);

parse(p, resGridSize, size, varargin{:});

% Handle constraints, to allow cell arrays as inputs
faceConstraints = p.Results.faceConstraints;
if iscell(faceConstraints)
    faceConstraints = constraintCellArrayToStruct(faceConstraints);
end
cellConstraints = p.Results.cellConstraints;
if iscell(cellConstraints)
    cellConstraints = constraintCellArrayToStruct(cellConstraints);
end
% Handle shape, to enable polygon shape
shape = p.Results.shape;
if length(shape) > 2
    shape = shapeArrayToStruct(shape);
end

py.gmsh4mrst.pebi_base_2D( ...
    cell_dimensions = p.Results.resGridSize, ...
    shape = shape, ...
    face_constraints = faceConstraints, ...
    face_constraint_factor = p.Results.faceConstraintFactor, ...
    face_constraint_parallel_factor = p.Results.faceConstraintParallelFactor, ...
    face_constraint_perpendicular_factor = p.Results.faceConstraintPerpendicularFactor, ...
    face_constraint_perpendicular_cells = p.Results.faceConstraintPerpendicularCells, ...
    face_constraint_point_factor = p.Results.faceConstraintPointFactor, ...
    face_constraint_refinement_factor = p.Results.faceConstraintRefinementFactor, ...
    min_FC_threshold_distance = p.Results.minFCThresholdDistance, ...
    max_FC_threshold_distance = p.Results.maxFCThresholdDistance, ...
    face_mesh_sampling = p.Results.faceMeshSampling, ...
    cell_constraints = cellConstraints, ...
    cell_constraint_factor = p.Results.cellConstraintFactor, ...
    cell_constraint_parallel_factor = p.Results.cellConstraintParallelFactor, ...
    cell_constraint_perpendicular_factor = p.Results.cellConstraintPerpendicularFactor, ...
    cell_constraint_perpendicular_cells = p.Results.cellConstraintPerpendicularCells, ...
    cell_constraint_point_factor = p.Results.cellConstraintPointFactor, ...
    cell_constraint_refinement_factor = p.Results.cellConstraintRefinementFactor, ...
    min_CC_threshold_distance = p.Results.minCCThresholdDistance, ...
    max_CC_threshold_distance = p.Results.maxCCThresholdDistance, ...
    CC_mesh_sampling = p.Results.CCMeshSampling, ...
    mesh_algorithm = p.Results.meshAlgorithm, ...
    recombination_algorithm = p.Results.recombinationAlgorithm, ...
    savename = "TEMP_Gmsh_MRST.m");

% Format shape again, this time for clippedPebi2D
bnd = p.Results.shape;
if length(bnd) == 2
    bnd = [0, 0; bnd(1), 0; bnd(1), bnd(2); 0, bnd(2)];
end

if isfile('TEMP_Gmsh_MRST.m')
    G = gmshToMRST('TEMP_Gmsh_MRST.m');
    if p.Results.convertToPEBI
        G = clippedPebi2D(G.nodes.coords, bnd);
    end
    delete TEMP_Gmsh_MRST.m
end

end


function constraintStruct = constraintCellArrayToStruct(constraints)
    if length(constraints) > 26
        error( ...
            "Constraints must have less than 26 separate lines to convert " ...
            + "automatically from cell array to Struct.")
    end
    constraintStruct = struct;
    for row = 1:length(constraints)
        rowChar = char(row + 96);    % +96 to start with char(a)
        constraintStruct.(rowChar).x = constraints{row}(:, 1);
        constraintStruct.(rowChar).y = constraints{row}(:, 2);
    end
end

function shapeStruct = shapeArrayToStruct(shape)
    shapeStruct = struct;
    shapeStruct.x = shape(:, 1);
    shapeStruct.y = shape(:, 2);
end