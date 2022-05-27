"""
Create a very simple GMSH plot, using user-supplied arguments, then save the
mesh to file.

Functions:
    pebi_grid_2D
"""
from array import array
from itertools import combinations
from math import atan2, sqrt, ceil
from typing import Any, Union, Iterable

import gmsh


def pebi_grid_2D(
        cell_dimensions: float,
        shape: list,
        *,
        face_constraints: Union[
                    'list[list[Iterable]]',
                    'dict[str, float]',
                    'dict[str, Iterable]',
                    'dict[str, dict[str, float]]',
                    'dict[Any, dict[str, Iterable]]'] = None,
        face_constraint_factor: float = 1/3,
        min_threshold_distance: float = 0.05,
        max_threshold_distance: float = 0.2,
        face_intersection_factor: float = None,
        min_intersection_distance: float = None,
        max_intersection_distance: float = None,
        fracture_mesh_sampling: int = 100,
        cell_constraints: Union[
                    'list[list[Iterable]]',
                    'dict[str, float]',
                    'dict[str, Iterable]',
                    'dict[str, dict[str, float]]',
                    'dict[Any, dict[str, Iterable]]'] = None,
        cell_constraint_factor: float = 1/4,
        cell_constraint_line_factor: float = None,
        cell_constraint_point_factor: float = None,
        mesh_algorithm: str = "Delaunay",
        recombination_algorithm: str = None,
        savename: str = "TEMP_Gmsh_MRST.m",
        run_frontend: bool = False
    ):
    """Create a 2D mesh, with user-supplied face constraints.

    This project was done as part of my Bachelor thesis during spring 2022, in
    order to help provide a secondary backend to the SINTEF-developed MATLAB
    module MRST, and its submodule UPR. For more information about MRST, see
    https://www.sintef.no/projectweb/mrst/. For more information about UPR, see
    https://www.sintef.no/projectweb/mrst/modules/upr/

    Args:
        cell_dimensions (float): Base dimensions of each cell.

        shape (Iterable[float] | Iterable[Iterable[float]] | dict[str, float]):
            Size or shape of the domain.
            Legal forms:
                Iterable[float]: Shape = (xmax, ymax), and is the size of a 
                    domain starting at (0, 0). Shape must have length 2.
                    Examples:
                        >>> shape = [1, 1]
                        >>> shape = (9001, 14)
                Iterable[Iterable[float]]: Shape is a list of points making up
                    the domain. All points in shape must be ordered either
                    clockwise or counterclockwise, and all points in shape must
                    have length >= 2 (only the first 2 elements are used).
                    Examples:
                        >>> shape = [
                            (0, 0), (0.5, 0.2), (1, 0), (1, 1), (0, 1)
                        ]
                        >>> shape = [
                            (0, 0), (1, 0), (1, 1), (0, 1)
                        ]
                dict[str, float]: Primarily used when complex domains are used
                    in MATLAB. dict must have keys "x" and "y", where the 
                    corresponding list being the x- and y-coordinates of the
                    complex domain.
                    Examples:
                        >>> shape = {
                                "x": [0, 0.5, 1, 1, 0],
                                "y": [0, 0.2, 0, 1, 1]
                            }
                        >>> shape = {"x": [0, 1, 1, 0], "y": [0, 0, 1, 1]}
            NOTE: Any constraints must be wholly within the supplied domain!

        face_constraints (list[Iterable] | dict[str, float] | dict[str, Iterable]
                | dict[str, dict[str, float]] | dict[str, dict[str, Iterable]],
                optional):
            Constraints for the cell faces. Each constraint is the coordinate
            of a surface-trace. The surface is assumed to be linear between 
            the coordinates. The resulting mesh will place sites such that 
            the surfaces are traced by faces of the grid.
            Legal forms:
                list[list[Iterable]]: Primarily used when calling the function
                    from Python. Each element of face constraint is a list of
                    point(s) along the constraint. Note that each element must
                    be a list of points, even if there is only one point!
                    Examples:
                        >>> face_constraints = [
                                [(0.1, 0.1), (0.9, 0.9)],
                                [(0.2, 0.8), (0.5, 0.5), (0.8, 0.2)]
                            ]
                        >>> face_constraints = [
                                [(0.5, 0.5)],
                                [(0.1, 0.1), (0.9, 0.9)]
                            ]
                dict[str, float] | dict[str, list]: Used for only one
                    constraint line, primarily when the function is called
                    from MATLAB. face_constraints must contain the keys "x" and
                    "y", with the corresponding float or list being the x- and
                    y-coordinates of point(s) along the constraint.
                    Examples:
                        >>> face_constraints = {"x": [0.1, 0.9], "y": [0.1, 0.9]}
                        >>> face_constraints = {"x": 0.5, "y": 0.5}
                dict[Any, dict[str, float]] | dict[Any, dict[str, Iterable]]:
                    Used for more than one constraint line, primarily when the
                    function is called from MATLAB. Each element in
                    face_constraints.values() must contain keys "x" and "y",
                    with the corresponding float or list being the x- and y-
                    coordinates of point(s) along the constraint. It does not
                    matter what is used as keys for the top-level dictionary.
                    Examples:
                        >>> face_constraints = {
                            "line": {"x": [0.1, 0.9], "y": [0.1, 0.9]},
                            "point": {"x": 0.5, "y": 0.5},
                            "keys dont matter": {"x": 0.1, "y": 0.9}.
                            2022: {"x": 0.9, "y": 0.1}
                        }
            NOTE: Any constraints must be wholly within the supplied domain!
            If None, face_constraints will be an empty list. Defaults to None.

        face_constraint_factor (float, optional): The size of the cells close
            to the face constraints, as compared to supplied cell_dimensions.
            Cells within min_threshold_distance will have size
            face_constraint_factor * cell_dimensions. Equivalent to FCFactor in
            MRST/UPR/pebiGrid2D. Defaults to 1/3.

        min_threshold_distance (float, optional): Distance from face constraints
            where cell dimensions will start increasing. Defaults to 0.05.

        max_threshold_distance (float, optional): Distance from face constraints
            where cell dimensions will be back to their default (max) value,
            i.e. the supplied argument cell_dimensions. Defaults to 0.2.

        face_intersection_factor (float, optional): The size of the cells close
            to intersections between face constraints, as compared to supplied
            cell_dimensions. Cells within min_intersection_distance from an 
            intersection will have size face_intersection_factor * cell_dimensions.
            If None, no extra cell shaping will occur around intersections.
            The factor is also used in "breaks" of lines, i.e. if there is a
            sharp "turn" in a line segment. Defaults to None.

        min_intersection_distance (float, optional): Distance from intersections
            where cell dimensions will start increasing. If None, will use
            min_threshold_distance. Defaults to None.

        max_intersection_distance (float, optional): Distance from intersections
            where cell dimensions will be back to their default (max) value,
            i.e. the supplied argument cell_dimensions. If None, will use
            min_threshold_distance. Defaults to None.

        fracture_mesh_sampling (int, optional): The number of points along the
            face constraints should be sampled to calculate the threshold
            distances. Defaults to 100.

        cell_constraints (list[Iterable] | dict[str, float] | dict[str, Iterable]
                | dict[str, dict[str, float]] | dict[str, dict[str, Iterable]],
                optional):
            Constraints for the cell centroids. Each constraint is a line,
            which the method will attempt to place in the center of the returned
            grid. The lines are assumed to be linear between the coordinates.
            If a constraint is only one coordinate, the line is treated as a
            point constraint.
            Legal forms:
                list[list[Iterable]]: Primarily used when calling the function
                    from Python. Each element of cell constraint is a list of
                    point(s) along the constraint. Note that each element must
                    be a list of points, even if there is only one point!
                    Examples:
                        >>> cell_constraints = [
                                [(0.1, 0.1), (0.9, 0.9)],
                                [(0.2, 0.8), (0.5, 0.5), (0.8, 0.2)]
                            ]
                        >>> cell_constraints = [
                                [(0.5, 0.5)],
                                [(0.1, 0.1), (0.9, 0.9)]
                            ]
                dict[str, float] | dict[str, list]: Used for only one
                    constraint line, primarily when the function is called
                    from MATLAB. cell_constraints must contain the keys "x" and
                    "y", with the corresponding float or list being the x- and
                    y-coordinates of point(s) along the constraint.
                    Examples:
                        >>> cell_constraints = {"x": [0.1, 0.9], "y": [0.1, 0.9]}
                        >>> cell_constraints = {"x": 0.5, "y": 0.5}
                dict[Any, dict[str, float]] | dict[Any, dict[str, Iterable]]:
                    Used for more than one constraint line, primarily when the
                    function is called from MATLAB. Each element in
                    cell_constraints.values() must contain keys "x" and "y",
                    with the corresponding float or list being the x- and y-
                    coordinates of point(s) along the constraint. It does not
                    matter what is used as keys for the top-level dictionary.
                    Examples:
                        >>> cell_constraints = {
                            "line": {"x": [0.1, 0.9], "y": [0.1, 0.9]},
                            "point": {"x": 0.5, "y": 0.5},
                            "keys dont matter": {"x": 0.1, "y": 0.9}.
                            2022: {"x": 0.9, "y": 0.1}
                        }
            NOTE: Any constraints must be wholly within the supplied domain!
            If None, cell_constraints will be an empty list. Defaults to None.

        cell_constraint_factor (float, optional): The size used for cells
            around cell constraints, as compared to supplied cell_dimensions.
            Equivalent to CCFactor in MRST/UPR/pebiGrid2D. Defaults to 1/4.

        cell_constraint_line_factor (float, optional): The size used for cells
            around cell constraint lines, as compared to the supplied
            cell_dimensions. Overrides cell_constraint_factor for lines. If
            set to None, cell_constraint_factor will be used for lines.
            Defaults to None.
            
        cell_constraint_point_factor (float, optional): The size used for cells
            around cell constraint points, as compared to the supplied
            cell_dimensions. Overrides cell_constraint_factor for points. If
            set to None, cell_constraint_factor will be used for points.
            Defaults to None.

        mesh_algorithm (str | int, optional): What meshing algorithm should be
            used. Can either be the Gmsh-given ID of the algorithm or a string:
                "MeshAdapt" = 1
                "Automatic" = 2
                "Delaunay" = 5
                "Frontal" = 6
                "BAMG" = 7
                "DelQuad" = 8
            Defaults to "Delaunay".

        recombination_algorithm (str | int, optional). What recombination
            algorithm should be used, and whether recombination should be done.
            Recombination makes Gmsh attempt to create a quadrangle mesh, rather
            than a triangle mesh. Can be either the Gmsh-given ID of the
            algorithm or a string:
                "Simple" = 0
                "Blossom" = 1
                "SimpleFull" = 2
                "BlossomFull" = 3
            If None, no recombination will be done.
            NOTE: Blossom is Gmsh default, but struggles with constraints. It
                may therefore be beneficial to use SimpleQuad, if there are
                constraints.
            NOTE 2: Recombination may lead to weird results for constraints,
                as the recombination is done after constraints have been applied.
                This is especially true for SimpleFull and BlossomFull, who
                automatically perform a coarser mesh, followed by recombination,
                smoothing and subdivision.
                
        savename (str, optional): Name of the saved file. If None, no file will
            be saved. The MATLAB functions assume that the file will be saved
            as "TEMP_Gmsh_MRST.m". Defaults to "TEMP_Gmsh_MRST.m".

        run_frontend (bool, optional): Set to True in order to run the Gmsh
            frontend and show the created mesh. Defaults to False.
    """
    if isinstance(shape, Iterable):
        if len(shape) < 2:
            raise ValueError(
                f"Shape must have length >= 2. Current length: {len(shape)}"
            )
        elif len(shape) == 2:
            # Shape is the size of the domain, starting at (0, 0)
            shape = [
                (0, 0), (0, shape[1]), (shape[0], shape[1]), (shape[0], 0)
            ]
    elif isinstance(shape, dict):
        # Shape is a dict, likely from MATLAB
        _assert_column_in_dict(shape, "x")
        _assert_column_in_dict(shape, "y")
        _assert_columns_have_same_length(shape, "x", "y")
        shape = [
            (x, y) for x, y in zip(shape.get("x"), shape.get("y"))
        ]
    else:
        raise ValueError(
            f"Type {type(shape)} is not supported for argument `shape`!"
        )

    if face_constraints is None:
        face_constraints = []
    if min_intersection_distance is None:
        min_intersection_distance = min_threshold_distance
    if max_intersection_distance is None:
        max_intersection_distance = max_threshold_distance
    if cell_constraints is None:
        cell_constraints = []
    if cell_constraint_line_factor is None:
        cell_constraint_line_factor = cell_constraint_factor
    if cell_constraint_point_factor is None:
        cell_constraint_point_factor = cell_constraint_factor

    # Handle mesh algorithm
    # According to 
    # https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t10.py
    # Frontal-Delaunay (6) usually leads to the highest quality meshes, but
    # Delaunay (5) handles complex mesh sizes better - especially size fields
    # with large element size gradients. We therefore default to Delaunay.
    # For quad-shaped grids, "Frontal-Delaunay for Quads" (8) may be beneficial.
    mesh_algorithm_dict = {
        "meshadapt": 1,
        "automatic": 2,
        "delaunay": 5,
        "frontal": 6,
        "bamg": 7,
        "delquad": 8
    }
    if isinstance(mesh_algorithm, int):
        if mesh_algorithm not in mesh_algorithm_dict.values():
            raise ValueError(
                "mesh_algorithm must be a legal value. Current value: "
                + f"{mesh_algorithm}. Legal values: {mesh_algorithm_dict}"
            )
    else:
        mesh_algorithm = mesh_algorithm_dict.get(mesh_algorithm.lower(), 5)
    
    # Do the same for the recombination algorithm
    recombination_algorithm_dict = {
        "simple": 0,
        "blossom": 1,
        "simplefull": 2,
        "blossomfull": 3
    }
    if isinstance(recombination_algorithm, int):
        if recombination_algorithm not in recombination_algorithm_dict.values():
            raise ValueError(
                "recombination_algorithm must be a legal value. Current value: "
                + f"{recombination_algorithm}. Legal values: {recombination_algorithm_dict}"
            )
    elif recombination_algorithm is None:
        pass
    else:
        recombination_algorithm = recombination_algorithm_dict.get(
            recombination_algorithm.lower()
        )
    
    # Do some (massive) argument handling, for use in MATLAB
    face_constraints = _format_constraints(face_constraints)
    cell_constraints = _format_constraints(cell_constraints)
    
    gmsh.initialize()
    gmsh.model.add("pebiGrid2D")

    # Create corners
    corners = [
        gmsh.model.geo.add_point(point[0], point[1], 0) for point in shape
    ]

    # Create circumference face_constraints
    circumference_face_constraints = [
        gmsh.model.geo.add_line(corners[i], corners[i+1])
            for i in range(len(corners) - 1)
    ] + [
        gmsh.model.geo.add_line(corners[-1], corners[0])
    ]

    # Create fractures (face constraints)
    fracture_points = []
    fractures = []
    if face_intersection_factor is not None:
        line_segments = []
    # Create fracture points
    for line in face_constraints:
        # Handle single points their own way
        if len(line) == 1:
            fracture_points.append(
                gmsh.model.geo.add_point(line[0][0], line[0][1], 0)
            )
            continue
        
        # There exists a line with at least 1 segment
        # As line segment [i] starts with the end of line segment [i - 1],
        # we avoid doubling on points by moving fracture_start out of the loop
        fracture_start = gmsh.model.geo.add_point(line[0][0], line[0][1], 0)
        for i in range(1, len(line)):
            fracture_end = gmsh.model.geo.add_point(line[i][0], line[i][1], 0)
            fractures.append(
                gmsh.model.geo.add_line(fracture_start, fracture_end)
            )
            fracture_start = fracture_end

            if face_intersection_factor is not None:
                line_segments.append([(line[i-1][0], line[i-1][1]), (line[i][0], line[i][1])])
    
    # Calculate intersections
    intersection_points = []
    if face_intersection_factor is not None:
        for line_1, line_2 in combinations(line_segments, 2):
            intersection = _find_intersection(
                line_1[0], line_1[1], line_2[0], line_2[1]
            )
            if intersection is not None:
                intersection_points.append(
                    gmsh.model.geo.add_point(intersection[0], intersection[1], 0)
                )

    # Create cell constraints
    cc_loops = []               # Holds line loops - used when making surface
    cc_point_surfaces = []      # Holds surfaces around CC points
    cc_line_surfaces = []       # Holds surfaces around CC lines
    cc_point_size = cell_constraint_point_factor * cell_dimensions  # Cell size around points
    cc_line_size = cell_constraint_line_factor * cell_dimensions    # Cell size around lines
    for line in cell_constraints:
        if len(line) == 1:
            # line is a single point
            x, y = line[0][0], line[0][1]
            # Create corners around the CC point - up, down, left and right
            surrounding_points = [
                gmsh.model.geo.add_point(x - cc_point_size/2, y, 0),
                gmsh.model.geo.add_point(x, y + cc_point_size/2, 0),
                gmsh.model.geo.add_point(x + cc_point_size/2, y, 0),
                gmsh.model.geo.add_point(x, y - cc_point_size/2, 0)
            ]
            # Create lines surrounding the CC point, i.e. between surrounding
            # points. These make up the mesh cell around the CC point
            surrounding_lines = [
                gmsh.model.geo.add_line(surrounding_points[0], surrounding_points[1]),
                gmsh.model.geo.add_line(surrounding_points[1], surrounding_points[2]),
                gmsh.model.geo.add_line(surrounding_points[2], surrounding_points[3]),
                gmsh.model.geo.add_line(surrounding_points[3], surrounding_points[0]),
            ]
            # Save the curve loop created
            cc_loops.append(
                gmsh.model.geo.add_curve_loop(surrounding_lines)
            )
            # Save the surface created
            cc_point_surfaces.append(gmsh.model.geo.add_plane_surface([cc_loops[-1]]))
            # Make each surrounding line into a transfinite curve with 2 points
            # This creates a point in each corner, leading to a single cell
            # within the surrounding lines - naturally with a face perfectly on
            # the CC point
            for sur_line in surrounding_lines:
                gmsh.model.geo.mesh.set_transfinite_curve(sur_line, 2)
            # Make the surface a transfinite surface
            gmsh.model.geo.mesh.set_transfinite_surface(cc_point_surfaces[-1])
            # Convert the (perfectly triangle) surface into a quadrangle one.
            gmsh.model.geo.mesh.set_recombine(2, cc_point_surfaces[-1])
        
        else:
            # line has at least 1 segment
            # We first handle the start, then all mid points, then the end
            # Compute the normal vector of the line from start- to next point
            delta_x = line[1][0] - line[0][0]
            delta_y = line[1][1] - line[0][1]
            normal_x, normal_y = _get_perpendicular(delta_x, delta_y)

            # Create starting points, one along the normal vector and one
            # opposite of the normal vector
            point_1, point_2 = _get_extruded_points(line[0], normal_x, normal_y, cc_line_size)
            
            # Create actual Gmsh points of the starting points
            start_1 = gmsh.model.geo.add_point(point_1[0], point_1[1], 0)
            start_2 = gmsh.model.geo.add_point(point_2[0], point_2[1], 0)
            # Make a start line for the polygon surrounding our CC line
            start_line = gmsh.model.geo.add_line(start_1, start_2)
            # Convert the start line into a transfinite curve. We again use 2
            # transfinite points, leading to a width of 1 cell
            gmsh.model.geo.mesh.set_transfinite_curve(start_line, 2)

            # Handle all midpoints
            for i in range(1, len(line) - 1):
                # We find the normal vector by finding the midpoint of point
                # [i-1] and [1+1], then the vector from this midpoint to [i]
                midpoint = (
                    0.5 * (line[i-1][0] + line[i+1][0]),
                    0.5 * (line[i-1][1] + line[i+1][1]),
                )
                normal_x = line[i][0] - midpoint[0]
                normal_y = line[i][1] - midpoint[1]
                # If the points are in a line, we simply use the normal vector
                # from earlier
                if normal_x == normal_y == 0:
                    normal_x, normal_y = normal_x, normal_y
                # Like for the start point, we create two extrusions - one 
                # along the normal vector, and one opposite of it
                end_point_1, end_point_2 = _get_extruded_points(
                    line[i], normal_x, normal_y, cc_line_size
                )
                # If the line bends towards the right (creating an A-shape),
                # then we must flip the points. This ensures that the start-
                # and end points make a nice rectangle (compared to the twisted
                # rectangle we'd get otherwise)
                if _line_bends_towards_right(line[i-1], line[i], line[i+1]):
                    end_point_1, end_point_2 = end_point_2, end_point_1
                
                # Create actual Gmsh points of the extrusion points
                end_1 = gmsh.model.geo.add_point(end_point_1[0], end_point_1[1], 0)
                end_2 = gmsh.model.geo.add_point(end_point_2[0], end_point_2[1], 0)
                
                # Compute the transfinite points of the parallel lines
                # We use the line from start_2 -> end_1 as our measuring stick, as
                # the difference between the parallel lines is likely very small
                distance = _distance(line[i-1], line[i])
                parallel_line_points = _calculate_parallel_transfinite_points(
                    distance, cc_line_size
                )

                end_line = _create_transfinite_cc_box(
                    start_1, start_2, start_line,
                    end_1, end_2,
                    cc_loops, cc_line_surfaces, parallel_line_points
                )

                # Convert the previous end line to a new start line
                # Note the flip (end_2 -> start_1 and vice versa)
                # This ensures we get neat rectangles for each line segment
                start_1, start_2 = end_2, end_1
                start_line = -end_line
            
            # Handle the end point
            # Same as for the start point - compute the normal vector
            delta_x = line[-1][0] - line[-2][0]
            delta_y = line[-1][1] - line[-2][1]
            normal_x, normal_y = _get_perpendicular(delta_x, delta_y)
            # Use the normal vector to get extrusion points
            end_point_2, end_point_1 = _get_extruded_points(line[-1], normal_x, normal_y, cc_line_size)

            # Create actual Gmsh points from the extrusion points
            end_1 = gmsh.model.geo.add_point(end_point_1[0], end_point_1[1], 0)
            end_2 = gmsh.model.geo.add_point(end_point_2[0], end_point_2[1], 0)

            # Compute the transfinite points of the parallel lines
            # We use the line from start_2 -> end_1 as our measuring stick, as
            # the difference between the parallel lines is likely very small
            distance = _distance(line[-1], line[-2])
            parallel_line_points = _calculate_parallel_transfinite_points(
                distance, cc_line_size
            )

            # Create a transfinite surface out of the created CC box
            _create_transfinite_cc_box(
                start_1, start_2, start_line,
                end_1, end_2,
                cc_loops, cc_line_surfaces, parallel_line_points
            )

    # Create curve loop of circumference
    circumference = gmsh.model.geo.add_curve_loop(circumference_face_constraints)
    
    # Define surface from circumference
    # We remove the loops created from cell constraints
    surface = gmsh.model.geo.add_plane_surface([circumference, *cc_loops])

    # Synchronize to prepare for embedding fracture face_constraints
    gmsh.model.geo.synchronize()

    # Embed all fractures
    gmsh.model.mesh.embed(1, fractures, 2, surface)
    gmsh.model.mesh.embed(0, fracture_points, 2, surface)

    # Refine mesh around fractures. We want a finer mesh around the fractures,
    # to ensure as close fit of faces as possible. In addition, this helps
    # handle intersecting fractures
    # This is gotten from
    # https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t10.py

    # The Distance field returns the distance to (sampling) points on each
    # fracture
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", fracture_points)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", fractures)
    gmsh.model.mesh.field.setNumber(1, "Sampling", fracture_mesh_sampling)

    # The Treshold field uses the value from the Distance field to define a
    # change in element size depending on the computed distances
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin",
        face_constraint_factor * cell_dimensions
    )
    gmsh.model.mesh.field.setNumber(2, "SizeMax", cell_dimensions)
    gmsh.model.mesh.field.setNumber(2, "DistMin", min_threshold_distance)
    gmsh.model.mesh.field.setNumber(2, "DistMax", max_threshold_distance)

    # Add field for intersection points
    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "PointsList", intersection_points)
    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", cell_dimensions)
    gmsh.model.mesh.field.setNumber(4, "DistMin", min_intersection_distance)
    gmsh.model.mesh.field.setNumber(4, "DistMax", max_intersection_distance)
    if face_intersection_factor is not None:
        gmsh.model.mesh.field.setNumber(4, "SizeMin",
            face_intersection_factor * cell_dimensions
        )

    # We use the minimum of all fields as our background mesh
    gmsh.model.mesh.field.add("Min", 10)
    gmsh.model.mesh.field.setNumbers(10, "FieldsList", [2, 4, ])
    gmsh.model.mesh.field.setAsBackgroundMesh(10)

    # As we define the entire element size in our background mesh, we disable
    # certain on-by-default mesh calculations
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # We finally set the meshing algorithm.
    gmsh.option.setNumber("Mesh.Algorithm", mesh_algorithm)

    if recombination_algorithm is not None:
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", recombination_algorithm)
        gmsh.model.mesh.set_recombine(2, surface)

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Optionally save
    if savename is not None:
        gmsh.write(savename)
    
    # Optionally run frontend
    if run_frontend:
        gmsh.fltk.run()

    # Always finalize
    gmsh.finalize()




############################## HELPER FUNCTIONS ##############################
def _assert_column_in_dict(face_constraints: dict, column: str):
    if column not in face_constraints.keys():
        raise ValueError(
            f"Dictionary {dict} must contain column `{column}`!"
        )

def _assert_columns_have_same_length(face_constraints: dict, col1: str, col2: str):
    if len(face_constraints.get(col1)) != len(face_constraints.get(col2)):
        raise ValueError(
            f"Column `{col1}` and `{col2}` must be the same length! "
            + f"len({col1}): {len(face_constraints.get(col1))}, len({col2}): {len(face_constraints.get(col2))}"
        )

def _check_array_dict_and_return_line_list(face_constraints: dict):
    # Check that face_constraints contains column x and y
    _assert_column_in_dict(face_constraints, "x")
    _assert_column_in_dict(face_constraints, "y")
    
    # Handle actual content
    if isinstance(face_constraints.get("x"), float):
        # If the 'line' is a single point, then x and y are floats,
        # and must be handled explicitly
        return [(face_constraints.get("x"), face_constraints.get("y"))]
    else:
        # Assume 'line' is a list/tuple/Iterable of floats
        _assert_columns_have_same_length(face_constraints, "x", "y")
        return list(zip(face_constraints.get("x"), face_constraints.get("y")))

def _find_intersection(line_1_start, line_1_end, line_2_start, line_2_end):
    """Find the intersection of two line segments.

    Method inspired by
    https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect#565282
    with credit, quote:
        "Credit: this method is the 2-dimensional specialization of the 3D line
        intersection algorithm from the article 'Intersection of two lines in
        three-space' by Ronald Goldman, published in Graphics Gems, page 304

    Args:
        line_1_start (tuple): Start point of line 1
        line_1_end (tuple): End point of line 1
        line_2_start (tuple): Start point of line 2
        line_2_end (tuple): End point of line 2
    """
    if line_1_start[0] > line_1_end[0]:
        line_1_start, line_1_end = line_1_end, line_1_start
    if line_2_start[0] > line_2_end[0]:
        line_2_start, line_2_end = line_2_end, line_2_start
    
    def cross_product(v, w):
        return v[0]*w[1] - v[1]*w[0]
    
    delta_1 = [line_1_end[0]-line_1_start[0], line_1_end[1]-line_1_start[1]] #r
    delta_2 = [line_2_end[0]-line_2_start[0], line_2_end[1]-line_2_start[1]] #s

    if cross_product(delta_1, delta_2) == 0:
        # Either parallel or collinear -> No intersections we need to consider
        return None
    
    line_difference = [                         # (q - p)
        line_2_start[0] - line_1_start[0],
        line_2_start[1] - line_1_start[1]
    ]
    t = cross_product(line_difference, delta_2) / cross_product(delta_1, delta_2)
    u = cross_product(line_difference, delta_1) / cross_product(delta_1, delta_2)
    if 0 <= t <= 1 and 0 <= u <= 1:
        return [line_1_start[0] + t * delta_1[0], line_1_start[1] + t * delta_1[1]]
    else:
        return None


def _format_constraints(constraints) -> list:
    if isinstance(constraints, dict):
        if len(constraints) == 0:
            # constraints is an empty list -> No fractures to handle,
            # simply create a single mesh within size
            constraints = []
        elif isinstance(list(constraints.values())[0], array):
            # If data is sent from MATLAB, then constraints is a dict
            # If there is only one line, then we accept it as a 1D dict, i.e.
            # one can write
            #   constraints.x = [x1 ...];
            #   constraints.y = [y1 ...];
            # Then constraints is a dict in the shape {x: array, y: array}
            constraints = [_check_array_dict_and_return_line_list(constraints)]

        elif isinstance(list(constraints.values())[0], dict):
            # If data is sent from MATLAB, then constraints is a dict of dicts
            # If there are more than one line, then it must be a 2D dict, i.e.
            # one can write
            #   constraints.line1.x = [x11 ...];
            #   constraints.line1.y = [y11 ...];
            #   constraints.line2.x = [x21 ...];
            #   constraints.line2.y = [y21 ...];
            # Each element in constraints.values is a dict of {x: array, y: array}
            new_constraints = []
            for constraint in constraints.values():
                new_constraints.append(
                    _check_array_dict_and_return_line_list(constraint)
                )
            constraints = new_constraints
    return constraints

def _get_perpendicular(delta_x, delta_y) -> 'tuple[float, float]':
    return -delta_y, delta_x

def _get_extruded_points(
            base_point, normal_x, normal_y, cc_size
        ) -> 'tuple[tuple[float, float], tuple[float, float]]':
    # "Normalize" normal_x and normal_y, such that the length of the normal
    # vector = 1
    prev_length = sqrt(normal_x**2 + normal_y**2)
    normal_x = normal_x / prev_length
    normal_y = normal_y / prev_length
    extruded_above = (
        base_point[0] + normal_x * cc_size / 4,
        base_point[1] + normal_y * cc_size / 4
    )
    extruded_below = (
        base_point[0] - normal_x * cc_size / 4,
        base_point[1] - normal_y * cc_size / 4
    )
    return extruded_above, extruded_below

def _line_bends_towards_right(start_point, mid_point, end_point) -> bool:
    angle_start_mid = atan2(
        mid_point[1] - start_point[1], mid_point[0] - start_point[0]
    )
    angle_start_end = atan2(
        end_point[1] - start_point[1], end_point[0] - start_point[0]
    )
    return angle_start_mid > angle_start_end

def _create_transfinite_cc_box(
        start_1,
        start_2,
        start_line,
        end_1,
        end_2,
        cc_loops: list,
        cc_line_surfaces: list,
        parallel_line_points: int
    ) -> int:
    """Create a transfinite box from the supplied points
    
    Add the curve loop to cc_loops, and the surface to cc_line_surfaces.

    Return end_line
    """
    # Create the parallel lines of the box surrounding the CC line
    parallel_line_1 = gmsh.model.geo.add_line(start_2, end_1)
    parallel_line_2 = gmsh.model.geo.add_line(end_2, start_1)
    # Create the end line
    end_line = gmsh.model.geo.add_line(end_1, end_2)

    # Save the curve loop, for use when creating the base surface
    curve_loop = gmsh.model.geo.add_curve_loop([
        start_line,
        parallel_line_1,
        end_line,
        parallel_line_2
    ])
    cc_loops.append(curve_loop)

    # And save the actual surface created
    surface = gmsh.model.geo.add_plane_surface([curve_loop])
    cc_line_surfaces.append(surface)

    # Create a transfinite surface out of the surface we made from
    # the box surrounding the CC line
    # End line always has 2 transfinite points, to get nice single-width cells
    gmsh.model.geo.mesh.set_transfinite_curve(end_line, 2)
    gmsh.model.geo.mesh.set_transfinite_curve(parallel_line_1, parallel_line_points)
    gmsh.model.geo.mesh.set_transfinite_curve(parallel_line_2, parallel_line_points)
    gmsh.model.geo.mesh.set_transfinite_surface(surface)
    gmsh.model.geo.mesh.set_recombine(2, surface)

    return end_line

def _distance(point_1, point_2) -> float:
    """Compute the distance between two points, using Euclidean distance"""
    return sqrt((point_2[1] - point_1[1])**2 + (point_2[0] - point_1[0])**2)

def _calculate_parallel_transfinite_points(distance, cc_line_size):
    return ceil(distance / cc_line_size) + 1

############################## END OF HELPERS ##############################



if __name__ == "__main__":
    pebi_grid_2D(
        0.2, 
        face_constraints=[
            [(0.25, 0.25), (0.4, 0.5), (0.7, 0.7)],
            [(0.8, 0.1), (0.9, 0.2)],
            [(0.2, 0.9), (0.9, 0.1)]
        ],
        cell_constraints=[
            [(0.1, 0.1)],
            [(0.5, 0.8), (0.6, 0.7), (0.7, 0.8), (0.9, 0.6)],
        ],
        shape=[
            (0, 0), (0.5, 0.2), (1, 0), (1, 1), (0, 1)
        ],
        face_constraint_factor = 1/3,
        face_intersection_factor = 1/9,
        mesh_algorithm="DelQuad",
        # recombination_algorithm="simplefull",
        savename=None,
        run_frontend=True)