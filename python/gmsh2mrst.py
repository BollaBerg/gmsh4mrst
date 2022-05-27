"""
Create a very simple GMSH plot, using user-supplied arguments, then save the
mesh to file.

Functions:
    pebi_grid_2D
"""
from itertools import combinations
from typing import Any, Union, Iterable

import gmsh

from _arguments import (
    format_constraints, format_shape,
    format_meshing_algorithm, format_recombination_algorithm
)
from _geometry import (
    find_intersection, get_perpendicular, get_extruded_points, get_midpoint,
    line_bends_towards_right, distance, calculate_number_of_points, 
)
from _gmsh import (
    create_transfinite_cc_box, create_threshold_field, create_circumference
)

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
    # Handle default values
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

    # Format Gmsh algorithms
    mesh_algorithm = format_meshing_algorithm(mesh_algorithm)
    recombination_algorithm = format_recombination_algorithm(recombination_algorithm)
    
    # Format all complex arguments
    shape = format_shape(shape)
    face_constraints = format_constraints(face_constraints)
    cell_constraints = format_constraints(cell_constraints)
    
    gmsh.initialize()
    gmsh.model.add("gmsh4mrst")

    # Create corners
    corners = [
        gmsh.model.geo.add_point(point[0], point[1], 0) for point in shape
    ]

    # Create circumference
    circumference = create_circumference(corners)

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
            intersection = find_intersection(
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
            surrounding_lines = create_circumference(surrounding_points)

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
            normal_x, normal_y = get_perpendicular(delta_x, delta_y)

            # Create starting points, one along the normal vector and one
            # opposite of the normal vector
            point_1, point_2 = get_extruded_points(line[0], normal_x, normal_y, cc_line_size)
            
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
                midpoint = get_midpoint(line[i-1], line[i+1])
                normal_x = line[i][0] - midpoint[0]
                normal_y = line[i][1] - midpoint[1]
                # If the points are in a line, we simply use the normal vector
                # from earlier
                if normal_x == normal_y == 0:
                    normal_x, normal_y = normal_x, normal_y
                # Like for the start point, we create two extrusions - one 
                # along the normal vector, and one opposite of it
                end_point_1, end_point_2 = get_extruded_points(
                    line[i], normal_x, normal_y, cc_line_size
                )
                # If the line bends towards the right (creating an A-shape),
                # then we must flip the points. This ensures that the start-
                # and end points make a nice rectangle (compared to the twisted
                # rectangle we'd get otherwise)
                if line_bends_towards_right(line[i-1], line[i], line[i+1]):
                    end_point_1, end_point_2 = end_point_2, end_point_1
                
                # Create actual Gmsh points of the extrusion points
                end_1 = gmsh.model.geo.add_point(end_point_1[0], end_point_1[1], 0)
                end_2 = gmsh.model.geo.add_point(end_point_2[0], end_point_2[1], 0)
                
                # Compute the transfinite points of the parallel lines
                # We use the line from start_2 -> end_1 as our measuring stick, as
                # the difference between the parallel lines is likely very small
                line_length = distance(line[i-1], line[i])
                parallel_line_points = calculate_number_of_points(
                    line_length, cc_line_size
                )

                end_line = create_transfinite_cc_box(
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
            normal_x, normal_y = get_perpendicular(delta_x, delta_y)
            # Use the normal vector to get extrusion points
            end_point_2, end_point_1 = get_extruded_points(line[-1], normal_x, normal_y, cc_line_size)

            # Create actual Gmsh points from the extrusion points
            end_1 = gmsh.model.geo.add_point(end_point_1[0], end_point_1[1], 0)
            end_2 = gmsh.model.geo.add_point(end_point_2[0], end_point_2[1], 0)

            # Compute the transfinite points of the parallel lines
            # We use the line from start_2 -> end_1 as our measuring stick, as
            # the difference between the parallel lines is likely very small
            line_length = distance(line[-1], line[-2])
            parallel_line_points = calculate_number_of_points(
                line_length, cc_line_size
            )

            # Create a transfinite surface out of the created CC box
            create_transfinite_cc_box(
                start_1, start_2, start_line,
                end_1, end_2,
                cc_loops, cc_line_surfaces, parallel_line_points
            )

    # Create curve loop of circumference
    circumference_loop = gmsh.model.geo.add_curve_loop(circumference)
    
    # Define surface from circumference
    # We remove the loops created from cell constraints
    surface = gmsh.model.geo.add_plane_surface([circumference_loop, *cc_loops])

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

    fracture_threshold = create_threshold_field(
        point_list=fracture_points,
        curve_list=fractures,
        sampling=fracture_mesh_sampling,
        min_size=face_constraint_factor * cell_dimensions,
        max_size=cell_dimensions,
        min_distance=min_threshold_distance,
        max_distance=max_threshold_distance
    )

    # Add field for intersection points
    if face_intersection_factor is not None:
        intersection_min_size = face_intersection_factor * cell_dimensions
    else:
        intersection_min_size = cell_dimensions
    intersection_threshold = create_threshold_field(
        point_list=intersection_points,
        curve_list=None,
        sampling=None,
        min_size=intersection_min_size,
        max_size=cell_dimensions,
        min_distance=min_intersection_distance,
        max_distance=max_intersection_distance
    )

    # We use the minimum of all fields as our background mesh
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [
        fracture_threshold, intersection_threshold
    ])
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

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