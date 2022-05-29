function [G, Pts, F] = pebiGrid2DGmsh(resGridSize, shape, varargin)
% Construct a 2D PEBI grid, using Gmsh.
%
% This method creates a background mesh using background_grid_2D, then
% creates constraints using UPR. It is a semi-direct copy of pebiGrid2D,
% but with Gmsh as a drop-in replacement of Distmesh.
% 
% SYNOPSIS:
%   G = pebiGrid2DGmsh(resGridSize)
%   G = pebiGrid2DGmsh(..., 'Name1', Value1, ...)
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
% OPTIONAL PARAMETERS   (gmsh4mrst)
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
%                   Defaults to 1.
%
%   minThresholdDistance - Float. Distance from face constraints where cell
%                   dimensions will start scaling size. Defaults to 0.05.
%
%   maxThresholdDistance - Float. Distance from face constraints where cell
%                   dimensions will be back to their default (max) value,
%                   i.e. the supplied argument resGridSize. Defaults
%                   to 0.2.
%
%   faceIntersectionFactor - Float. The size of the cells close to
%                   intersections between face constraints, as compared to
%                   resGridSize. Cells within minIntersectionDistance
%                   from an intersection will have size
%                   faceIntersectionFactor * resGridSize. If missing, no
%                   extra cell shaping will occur around intersections.
%                   Defaults to missing (aka python None).
%
%   minIntersectionDistance - Float. Distance from intersections where cell
%                   dimensions will start increasing. If missing, 
%                   minThresholdDistance will be used. Defaults to missing.
%
%   maxIntersectionDistance - Float. Distance from intersections where
%                   cell dimensions will be back to their default (max)
%                   value, i.e. the supplied argument resGridSize. If
%                   missing, minThresholdDistance will be used. Defaults to
%                   missing.
%
%   fractureMeshSampling - Int. The number of points along the face
%                   constraints should be sampled to calculate the
%                   threshold distances. Defaults to 100.
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
%                   Defaults to 1.
%
%   cellConstraintLineFactor - Float. The size used for cells around cell
%                   constraint lines, as compared to the supplied
%                   resGridSize. Overrides cellConstraintFactor for
%                   lines. If missing, cellConstraintFactor will be
%                   used for lines. Defaults to missing.
%
%   cellConstraintPointFactor - Float. The size used for cells around cell
%                   constraint points, as compared to the supplied 
%                   resGridSize. Overrides cellConstraintFactor for
%                   points. If missing, cellConstraintFactor will be used
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
% OPTIONAL PARAMETERS   (pebiGrid2D)
% These are directly from pebiGrid2D, and do the same thing as there.
%
%   interpolateCC - logical. If false, each segment in the corresponding
%                   constraint line will be represented by at least one
%                   cell. If true, the routine will interpolate along the
%                   curve, which means that sites will not
%                   necessarily fall exactly on the prescribed curve. Can
%                   be either one element, or an array of one element per
%                   cell constraint. Defaults to false.
%
%   CCFactor - float. Relative grid size of the cell constraint grid cells
%                   compared to supplied resGridSize. If CCFactor = 0.5,
%                   then the constrained cells will be about half the size
%                   of the background cells. Defaults to 1/4.
%
%   CCRefinement - logical. Whether refinement should be done around the
%                   cell constraints. Defaults to false.
%
%   CCRho - Function. Gives the relative distance between cell constraint
%                   sites. If CCRho = 0.5 in an area the distance between
%                   constraint sites will be around half the default, i.e.
%                   around 0.5 * CCFactor * resGridSize. Defaults to
%                   @(x) ones(size(x,1),1)
%
%   protLayer - logical. Whether a protection layer is added on both sides
%                   of the cell constraints. Defaults to true.
%
%   protD - cell array of functions. Gives the distance from the cell
%                   constrained sites to the protection sites. The function
%                   is evaluated along the well path such that protD(0) is
%                   the start of the well while protD(1) is the end of the
%                   well. The array should have either one function or one
%                   function for each well path. Defaults to
%                   {@(x) ones(size(x,1),1)*resGridSize/4}
%
%   interpolateFC - logical. If false, each segment in the corresponding
%                   constraint line will be represented by at least one
%                   cell. If true, the routine will interpolate along the
%                   curve, which means that faces will not
%                   necessarily fall exactly on the prescribed curve. Can
%                   be either one element, or an array of one element per
%                   face constraint. Defaults to false.
%
%   FCFactor - float. Relative grid size of the surface grid cells compared
%                   to background cells. If FCFactor = 0.5, then the
%                   surface cells will be about half the size of the
%                   reservoir cells. Defaults to 1/2.
%
%   circleFactor - float. The ratio between the radius and distance between
%                   the circles used to create the face constraints. A
%                   small value will place the surface points close to the
%                   surfaces, while a large value will place them far from
%                   the surfaces. Valid values are between 0.5 and 1.
%                   Defaults to 0.6.
%
%   FCRho - Function. Gives the relative distance between
%                   surface sites along a face constraint. If FCRho = 0.5
%                   in an area, the face constraint sites will be about 50%
%                   closer than in other areas. Defaults to
%                   @(x) ones(size(x,1),1).
%
%   sufFCCond - logical. Whether we should enforce the sufficient and
%                   necessary face constraint condition. If false, we
%                   instead enforce a less strict condition and remove any
%                   background sites that are closer to the constraint
%                   sites than the constraint grid size. The conformity is
%                   then not guaranteed. May be set to false to handle
%                   problems with bad cells at the end of faults due to the
%                   sufficient condition removing some reservoir points.
%                   Defaults to true.
%
% RETURNS
%   G       - Valid grid definition.
%               The fields
%                   - G.cells.tag is TRUE for all cell constraints.
%                   - G.faces.tag is TRUE for all face constraints.
%   Pts     - Array [G.cells.num x 3] of the Pebi sites.
%   F       - Struct with elements as returned from surfaceSites2D
%
% EXAMPLE:
%   resGridSize = 0.1;
%   domain = [0 0; 0.5 0.2; 1 0; 1 1; 0 1;];
%   faceConstraints = {
%       [0.25 0.25; 0.4 0.5; 0.7 0.7;], ...
%       [0.8 0.1; 0.9 0.2], ...
%       [0.2 0.9; 0.9 0.1], ...
%   };
%   cellConstraints = {
%       [0.1 0.1], ...
%       [0.5 0.8; 0.6 0.7; 0.7 0.8; 0.9 0.6];
%   };
%   G = pebiGrid2DGmshBackground( ...
%       resGridSize, ...
%       domain, ...
%       'faceConstraints', faceConstraints, ...
%       'cellConstraints', cellConstraints, ...
%       'protLayer', true ...
%   );
%   plotGrid(G, 'faceColor', 'none');
%

defaultFaceConstraints = {};
defaultFaceConstraintFactor = 1;
defaultMinThresholdDistance = 0.05;
defaultMaxThresholdDistance = 0.2;
defaultFaceIntersectionFactor = string(missing);   % => Python None
defaultMinIntersectionDistance = string(missing);
defaultMaxIntersectionDistance = string(missing);
defaultFractureMeshSampling = 100;
defaultCellConstraints = {};
defaultCellConstraintFactor = 1;
defaultCellConstraintLineFactor = string(missing);
defaultCellConstraintPointFactor = string(missing);
defaultMeshAlgorithm = "Delaunay";
defaultRecombinationAlgorithm = string(missing);

% pebiGrid2D arguments
defaultInterpolateCC = false;
defaultCCFactor = 1/4;
defaultCCRefinement = false;
defaultCCRho = @(x) ones(size(x,1),1);
defaultProtLayer = false;
defaultProtD = {@(x) ones(size(x,1),1)*resGridSize/4};
defaultInterpolateFC = false;
defaultFCFactor = 1/2;
defaultCircleFactor = 0.6;
defaultFCRho = @(x) ones(size(x,1),1);
defaultSufFCCond = true;

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
    validShape = @(x) isnumeric(x) && length(x) >= 2;
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
addRequired(p, 'shape', validShape);
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

% pebiGrid2D parameters
addParameter(p, 'interpolateCC', defaultInterpolateCC, @islogical);
addParameter(p, 'CCFactor', defaultCCFactor, validFloat);
addParameter(p, 'CCRefinement', defaultCCRefinement, @islogical);
addParameter(p, 'CCRho', defaultCCRho);
addParameter(p, 'protLayer', defaultProtLayer, @islogical);
addParameter(p, 'protD', defaultProtD);
addParameter(p, 'interpolateFC', defaultInterpolateFC, @islogical);
addParameter(p, 'FCFactor', defaultFCFactor, validFloat);
addParameter(p, 'circleFactor', defaultCircleFactor, validFloat);
addParameter(p, 'FCRho', defaultFCRho);
addParameter(p, 'sufFCCond', defaultSufFCCond, @islogical);

parse(p, resGridSize, shape, varargin{:});
params = p.Results;

% Handle constraints, to allow cell arrays as inputs
faceConstraints = params.faceConstraints;
if iscell(faceConstraints)
    pyFaceConstraints = constraintCellArrayToStruct(faceConstraints);
end
cellConstraints = params.cellConstraints;
if iscell(cellConstraints)
    pyCellConstraints = constraintCellArrayToStruct(cellConstraints);
end

% Handle shape, to enable polygon shape
pyShape = params.shape;
if length(pyShape) > 2
    pyShape = shapeArrayToStruct(shape);
end

% Call Python, to compute background grid
py.gmsh4mrst.background_grid_2D(...
    cell_dimensions = params.resGridSize, ...
    shape = pyShape, ...
    face_constraints = pyFaceConstraints, ...
    face_constraint_factor = params.faceConstraintFactor, ...
    min_FC_threshold_distance = params.minThresholdDistance, ...
    max_FC_threshold_distance = params.maxThresholdDistance, ...
    face_intersection_factor = params.faceIntersectionFactor, ...
    min_intersection_distance = params.minIntersectionDistance, ...
    max_intersection_distance = params.maxIntersectionDistance, ...
    fracture_mesh_sampling = params.fractureMeshSampling, ...
    cell_constraints = pyCellConstraints, ...
    cell_constraint_factor = params.cellConstraintFactor, ...
    cell_constraint_line_factor = params.cellConstraintLineFactor, ...
    cell_constraint_point_factor = params.cellConstraintPointFactor, ...
    min_CC_threshold_distance = params.minThresholdDistance, ...
    max_CC_threshold_distance = params.maxThresholdDistance, ...
    mesh_algorithm = params.meshAlgorithm, ...
    recombination_algorithm = params.recombinationAlgorithm, ...
    savename = "TEMP_Gmsh_MRST.m");

if isfile("TEMP_Gmsh_MRST.m")
    G = gmshToMRST("TEMP_Gmsh_MRST.m");
    delete TEMP_Gmsh_MRST.m
else
    error("No file generated. Something failed in Python!")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% START OF CODE ADAPTED FROM pebiGrid2D %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adopt some variables used below
FCGridSize = params.FCFactor * params.resGridSize;
CCGridSize = params.FCFactor * params.resGridSize;
CCRef = params.CCRefinement;
circleFactor = params.circleFactor;
CCRho = @(x) CCGridSize * params.CCRho(x);
FCRho = params.FCRho;

% Set domain function
shape = params.shape;
if length(shape) == 2
    shape = [0, 0; shape(1), 0; shape(1), shape(2); 0, shape(2)];
end

% Format interpolateCC
if ~isempty(params.cellConstraints)
    if (numel(params.interpolateCC) == 1)
        % Extend interpolateCC to be an array of length = len(cellConstraints)
        params.interpolateCC = repmat(params.interpolateCC, numel(cellConstraints),1);
    end
    % Assert that each cellConstraint has an interpolateCC-value
    assert(numel(params.interpolateCC)==numel(cellConstraints));
    
    if numel(params.protD) == 1
        % Extend protD to be an array of length = len(cellConstraints)
        params.protD = repmat(params.protD,numel(cellConstraints),1);
    end
    % Assert that each cellConstraint has a protD
    assert(numel(params.protD) == numel(cellConstraints));
end

% Format interpolateFC
if ~isempty(params.faceConstraints)
    if (numel(params.interpolateFC) == 1)
        % Extend interpolateFC to be an array of length = len(faceConstraints)
        params.interpolateFC = repmat(params.interpolateFC, numel(params.faceConstraints),1);
    end
    % Assert that each faceConstraint has an interpolateFC
    assert(numel(params.interpolateFC)==numel(params.faceConstraints));
end

% Split face constraints
[faceConstraints, fCut, fcCut, IC] = splitAtInt2D(params.faceConstraints, params.cellConstraints);
interpFL = params.interpolateFC(IC);
% Split cell constraints
[cellConstraints,  cCut, cfCut, IC] = splitAtInt2D(params.cellConstraints, params.faceConstraints);
interpWP = params.interpolateCC(IC);
protD    = params.protD(IC);

% Find vertical cell constraints
nw              = cellfun(@numel, params.cellConstraints)/2;
vW              = nw==1;
cellConstraints = [cellConstraints,params.cellConstraints(vW)];
cCut            = [cCut;zeros(sum(vW),1)];
cfCut           = [cfCut; zeros(sum(vW),1)];
interpWP        = [interpWP; params.interpolateCC(vW)];
protD           = [protD; params.protD(vW)];


% Create cell constraint sites (wells)
sePtn = [cfCut==2|cfCut==3, cfCut==1|cfCut==3];
if CCRef
  fLen = CCGridSize*params.faceConstraintFactor*1.2;  
else
  fLen = FCGridSize;
end
bisectPnt = (fLen.^2 - (circleFactor*fLen).^2 ...
            + (circleFactor*fLen).^2) ./(2*fLen);
faultOffset = sqrt((circleFactor*fLen).^2 - bisectPnt.^2);
sePtn = (1.0+faultOffset/CCGridSize)*sePtn;

[wellPts, ~, protectionPts, ~] ...
    = lineSites2D(cellConstraints, CCGridSize, ...
                         'sePtn',         sePtn,...
                         'cCut',          cCut, ...
                         'protLayer',     params.protLayer,...
                         'protD',         protD, ...
                         'CCRho',         CCRho, ...
                         'interpolateCC', interpWP);

% Create distance functions
if CCRef && ~isempty(wellPts)
    hresw  = @(x) min((ones(size(x,1),1)/CCFactor), ...
        1.2*exp(minPdist2(x,wellPts)/CCEps));
    hfault = @(x) CCGridSize*params.faceConstraintFactor*hresw(x).*FCRho(x);
else
    hfault = @(p) FCGridSize*FCRho(p);
end

% Create surface sites (face constraints, faults)
F = surfaceSites2D(faceConstraints, FCGridSize,...
                          'circleFactor',  circleFactor,...
                          'fCut',          fCut, ...
                          'fcCut',         fcCut, ...
                          'interpolateFC', interpFL, ...
                          'distFun',       hfault);

% Remove tip sites outside domain
if size(F.t.pts, 1) > 0
    innside = inpolygon(F.t.pts(:, 1), F.t.pts(:, 2), shape(:, 1), shape(:, 2));
    F.t.pts = F.t.pts(innside, :);
end 

% Remove conflict points
if params.sufFCCond
	Pts = surfaceSufCond2D(G.nodes.coords,F);
else
	Pts = removeConflictPoints(G.nodes.coords,F.f.pts, F.f.Gs);
end
Pts  = [F.f.pts; wellPts; protectionPts; Pts];
% For some reason, we sometimes get imaginary points (although with the
% imaginary element = 0). We only care about real points!
% We also remove some possible duplicates
Pts = unique(real(Pts), 'rows');


% Create grid
G = clippedPebi2D(Pts, shape);

% label face constraint faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1; 
  % N == 1 is now a boundary constraint, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping for no surface pts.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num,1);
end

% Label cell constraint cells
if ~isempty(wellPts)
  G.cells.tag = false(G.cells.num,1);
  % Add tag to all cells generated from wellPts
  wellCells = size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1);
  G.cells.tag(wellCells)= true;
  
  % Add tag to cell and face constraint crossings
  endOfLine = fcCut==1 | fcCut==3;        % Crossing at end of face constraint
  strOfLine = fcCut==2 | fcCut==3;        % Crossing at start of face constraint
  fIde      = F.l.fPos([false;endOfLine]);
  fIds      = F.l.fPos(strOfLine);
  fToTag    = [F.l.f(mcolon(fIde - 2,fIde-1)); F.l.f(mcolon(fIds,fIds+1))];

  G.cells.tag(fToTag) = true;
else
  G.cells.tag = false(G.cells.num,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  END OF CODE ADAPTED FROM pebiGrid2D  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
