function [G, Pts, F] = pebiGrid2DGmshBackground(resGridSize, shape, varargin)
% Construct a 2D PEBI grid, using Gmsh.
%
% This method creates a background mesh using background_grid_2D, then
% creates constraints using UPR. It is a semi-direct copy of pebiGrid2D,
% but with Gmsh as a drop-in replacement of Distmesh.
% 
% SYNOPSIS:
%   G = pebiGrid2DGmsh2(resGridSize, varargin)
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
defaultProtD = {@(x) ones(size(x,1),1)*resGridSize/4};
defaultProtLayer = false;
defaultInterpolateFC = false;
defaultCCRefinement = false;
defaultCircleFactor = 0.6;
defaultFCRefinement = false;
defaultSufFCCond = true;
defaultUseMrstPebi = false;
defaultFCFactor = 1/2;
defaultCCFactor = 1/4;


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
addParameter(p, 'protD', defaultProtD);
addParameter(p, 'protLayer', defaultProtLayer, @islogical);
addParameter(p, 'interpolateFC', defaultInterpolateFC, @islogical);
addParameter(p, 'CCRefinement', defaultCCRefinement, @islogical);
addParameter(p, 'circleFactor', defaultCircleFactor, validFloat);
addParameter(p, 'FCRefinement', defaultFCRefinement, @islogical);
addParameter(p, 'sufFCCond', defaultSufFCCond, @islogical);
addParameter(p, 'useMrstPebi', defaultUseMrstPebi, @islogical);
addParameter(p, 'FCFactor', defaultFCFactor, validFloat);
addParameter(p, 'CCFactor', defaultCCFactor, validFloat);

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
shape = params.shape;
if length(shape) > 2
    pyShape = shapeArrayToStruct(shape);
end

% Set domain function
if length(shape) == 2
	rectangle = [0,0; shape(1), shape(2)];
	fd = @(p,varargin) drectangle(p, 0, shape(1), 0, shape(2));
	corners = [0,0; 0,shape(2); shape(1),0; shape(1), shape(2)];
	vararg  = [];
    shape = [0, 0; shape(1), 0; shape(1), shape(2); 0, shape(2)];
else
	rectangle = [min(shape); max(shape)];
	corners   = shape;
	fd        = @dpoly;
	vararg    = [shape; shape(1,:)];
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
FCRef = params.FCRefinement;
% Workaound on Rho functions, to get things working
CCRhoFunc = @(x) ones(size(x,1),1);
CCRho = @(x) CCGridSize*CCRhoFunc(x);
FCRho = CCRhoFunc;

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
    hresw  = @(p) constFunc(p)/CCFactor;
    hfault = @(p) FCGridSize*FCRho(p);
end

% Create surface sites (face constraints, faults)
F = surfaceSites2D(faceConstraints, FCGridSize,...
                          'circleFactor',  circleFactor,...
                          'fCut',          fCut, ...
                          'fcCut',         fcCut, ...
                          'interpolateFC', interpFL, ...
                          'distFun',       hfault);

if FCRef && ~isempty(F.f.pts)
  hresf = @(x) min((ones(size(x,1),1)/params.faceConstraintFactor), ...
                    1.2*exp(minPdist2(x,F.f.pts)/FCEps));
else
  hresf = @(p) constFunc(p)/CCFactor;
end

% Remove tip sites outside domain
if size(F.t.pts, 1) > 0
    innside = inpolygon(F.t.pts(:, 1), F.t.pts(:, 2), shape(:, 1), shape(:, 2));
    F.t.pts = F.t.pts(innside, :);
end 

if FCRef && CCRef
    ds   = min(CCGridSize,FCGridSize);
    hres = @(x,varargin) min(hresf(p), hresw(p));
elseif FCRef
    ds   = FCGridSize;
    hres = @(p,varargin) hresf(p);
else 
    ds   = CCGridSize;
    hres = @(p, varargin) hresw(p);
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
if params.useMrstPebi
	t    = pyOutput(Pts);
	% Fix boundary
	pmid = (Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;% Compute centroids
	t    = t(fd(pmid,vararg)<-0.001*CCFactor,:);   % Keep interior triangles

	G = triangleGrid(Pts, t);
	G = pebi(G);
else
    G = clippedPebi2D(Pts, shape);
end

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
