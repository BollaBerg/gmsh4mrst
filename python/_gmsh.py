import gmsh

def create_transfinite_cc_box(
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

def create_threshold_field(
        point_list: list,
        curve_list: list,
        sampling: int,
        min_size: float,
        max_size: float,
        min_distance: float,
        max_distance: float) -> int:
    """Create a Threshold field, using a Distance field as input

    Args:
        point_list (list): List of point IDs to calculate threshold for
        curve_list (list): List of curve IDs to calculate threshold for
        sampling (int): Number of sample points along the curves
        min_size (float): Minimum size of cells (close to points/curves)
        max_size (float): Maximum size of cells (far from points/curves)
        min_distance (float): Distance where ramp-up should start
        max_distance (float): Distance where ramp-up is finished

    Returns:
        int: ID of Threshold field
    """
    
    # The Distance field returns the distance to (sampling) points on each
    # fracture
    distance = gmsh.model.mesh.field.add("Distance")

    if point_list is not None:
        gmsh.model.mesh.field.setNumbers(distance, "PointsList", point_list)
    if curve_list is not None:
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", curve_list)
        gmsh.model.mesh.field.setNumber(distance, "Sampling", sampling)

    # The Treshold field uses the value from the Distance field to define a
    # change in element size depending on the computed distances
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", min_size)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", max_size)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", min_distance)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", max_distance)

    return threshold

def create_circumference(points: 'list[int]') -> 'list[int]':
    """Create a Gmsh circumference around a list of points.

    Args:
        points (list[int]): List of point ID to encapsulate

    Returns:
        list[int]: List of Gmsh IDs of the resulting lines
    """
    output = []
    for i in range(len(points) - 1):
        output.append(gmsh.model.geo.add_line(points[i], points[i+1]))
    output.append(gmsh.model.geo.add_line(points[-1], points[0]))
    return output

def create_fracture_point(point: tuple, intersection_IDs: dict) -> int:
    intersection_ID = intersection_IDs.get(point, "DEFAULT")
    if intersection_ID == "DEFAULT":
        return gmsh.model.geo.add_point(point[0], point[1], 0)
    elif intersection_ID is None:
        gmsh_point = gmsh.model.geo.add_point(point[0], point[1], 0)
        intersection_IDs[point] = gmsh_point
        return gmsh_point
    else:
        return intersection_IDs.get(point)

def create_cell_constraint_point(point: tuple,
                                cc_point_size: float,
                                cc_loops: list,
                                cc_point_surfaces: list):
    """Handle everything needed and create transfinite "box" for CC point

    Args:
        point (tuple): The CC point, with shape (x, y)
        cc_point_size (float): The diameter of the transfinite box
        cc_loops (list): A list of CC curve loops, to append to
        cc_point_surfaces (list): A list of CC point surfaces, to append to
    """
    x, y = point[0], point[1]
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