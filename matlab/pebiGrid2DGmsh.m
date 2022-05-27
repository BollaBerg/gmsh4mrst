function G = pebiGrid2DGmsh(resGridSize, size, varargin)
% Construct a 2D grid, using Gmsh
% 
% SYNOPSIS:
%   G = pebiGrid2DGmsh(resGridSize, pdims, faceConstraints)
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

defaultFaceConstraints = struct;
defaultFaceConstraintFactor = 1/3;
defaultMinThresholdDistance = 0.05;
defaultMaxThresholdDistance = 0.2;
defaultFaceIntersectionFactor = string(missing);   % => Python None
defaultMinIntersectionDistance = string(missing);
defaultMaxIntersectionDistance = string(missing);
defaultFractureMeshSampling = 100;
defaultCellConstraints = struct;
defaultCellConstraintFactor = 1/4;
defaultCellConstraintLineFactor = string(missing);
defaultCellConstraintPointFactor = string(missing);
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
addRequired(p, 'size', validSize);
addParameter(p, 'faceConstraints', defaultFaceConstraints);
addParameter(p, 'faceConstraintFactor', defaultFaceConstraintFactor, validFloat);
addParameter(p, 'minThresholdDistance', defaultMinThresholdDistance, validFloat);
addParameter(p, 'maxThresholdDistance', defaultMaxThresholdDistance, validFloat);
addParameter(p, 'faceIntersectionFactor', defaultFaceIntersectionFactor, @validOptionalFloat);
addParameter(p, 'minIntersectionDistance', defaultMinIntersectionDistance, @validOptionalFloat);
addParameter(p, 'maxIntersectionDistance', defaultMaxIntersectionDistance, @validOptionalFloat);
addParameter(p, 'fractureMeshSampling', defaultFractureMeshSampling, validInt);
addParameter(p, 'cellConstraints', defaultCellConstraints);
addParameter(p, 'cellConstraintFactor', defaultCellConstraintFactor, validFloat);
addParameter(p, 'cellConstraintLineFactor', defaultCellConstraintLineFactor, @validOptionalFloat);
addParameter(p, 'cellConstraintPointFactor', defaultCellConstraintPointFactor, @validOptionalFloat);
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
shape = p.Results.size;
if length(shape) > 2
    shape = shapeArrayToStruct(shape);
end

py.gmsh2mrst.pebi_grid_2D( ...
    cell_dimensions = p.Results.resGridSize, ...
    shape = shape, ...
    face_constraints = faceConstraints, ...
    face_constraint_factor = p.Results.faceConstraintFactor, ...
    min_threshold_distance = p.Results.minThresholdDistance, ...
    max_threshold_distance = p.Results.maxThresholdDistance, ...
    face_intersection_factor = p.Results.faceIntersectionFactor, ...
    min_intersection_distance = p.Results.minIntersectionDistance, ...
    max_intersection_distance = p.Results.maxIntersectionDistance, ...
    fracture_mesh_sampling = p.Results.fractureMeshSampling, ...
    cell_constraints = cellConstraints, ...
    cell_constraint_factor = p.Results.cellConstraintFactor, ...
    cell_constraint_line_factor = p.Results.cellConstraintLineFactor, ...
    cell_constraint_point_factor = p.Results.cellConstraintPointFactor, ...
    mesh_algorithm = p.Results.meshAlgorithm, ...
    recombination_algorithm = p.Results.recombinationAlgorithm, ...
    savename = "TEMP_Gmsh_MRST.m");
G = 0;

if isfile('TEMP_Gmsh_MRST.m')
    G = gmshToMRST('TEMP_Gmsh_MRST.m');
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