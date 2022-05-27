"""
    Helper files for internal use - Gmsh related.
    Copyright (C) 2022  Andreas Bjelland Berg

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import gmsh

from gmsh4mrst._geometry import (
    get_perpendicular, get_extruded_points, get_midpoint,
    line_bends_towards_right, distance, calculate_number_of_points,
)

def create_transfinite_cc_box(
        start_1,
        start_2,
        start_line,
        end_1,
        end_2,
        cc_loops: list,
        cc_line_surfaces: list,
        parallel_line_points: int,
        cc_lines: list,
        perp_nodes: int = 2
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

    cc_lines.extend([parallel_line_1, parallel_line_2, end_line])

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
    gmsh.model.geo.mesh.set_transfinite_curve(end_line, perp_nodes)
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
                                cc_point_surfaces: list,
                                cc_lines: list,
                                nodes: int = 2):
    """Handle everything needed and create transfinite "box" for CC point

    Args:
        point (tuple): The CC point, with shape (x, y)
        cc_point_size (float): The diameter of the transfinite box
        cc_loops (list): A list of CC curve loops, to append to
        cc_point_surfaces (list): A list of CC point surfaces, to append to
        cc_lines (list): A list of CC lines, to append to
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
    cc_lines.extend(surrounding_lines)

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
        gmsh.model.geo.mesh.set_transfinite_curve(sur_line, nodes)
    # Make the surface a transfinite surface
    gmsh.model.geo.mesh.set_transfinite_surface(cc_point_surfaces[-1])
    # Convert the (perfectly triangle) surface into a quadrangle one.
    gmsh.model.geo.mesh.set_recombine(2, cc_point_surfaces[-1])


def create_cell_constraint_line(line: 'list[tuple]',
                                cc_parallel_size: float,
                                cc_perpendicular_size: float,
                                cc_loops: list,
                                cc_line_surfaces: list,
                                cc_lines: list,
                                perp_nodes: int = 2):
    """Handle everything needed and create transfinite "box" for CC line


    Args:
        line (list[tuple]): The CC line, in shape [(x1, y1), (x2, y2), ...]
        cc_line_size (float): The wanted length of each line segment
        cc_loops (list): A list of CC curve loops, to append to
        cc_line_surfaces (list): A list of CC line surfaces, to append to
        cc_lines (list): A list of CC lines, to append to
    """
    # We first handle the start, then all mid points, then the end
    # Compute the normal vector of the line from start- to next point
    delta_x = line[1][0] - line[0][0]
    delta_y = line[1][1] - line[0][1]
    normal_x, normal_y = get_perpendicular(delta_x, delta_y)

    # Create starting points, one along the normal vector and one
    # opposite of the normal vector
    point_1, point_2 = get_extruded_points(line[0], normal_x, normal_y, cc_perpendicular_size)
    
    # Create actual Gmsh points of the starting points
    start_1 = gmsh.model.geo.add_point(point_1[0], point_1[1], 0)
    start_2 = gmsh.model.geo.add_point(point_2[0], point_2[1], 0)
    # Make a start line for the polygon surrounding our CC line
    start_line = gmsh.model.geo.add_line(start_1, start_2)
    cc_lines.append(start_line)
    # Convert the start line into a transfinite curve. We again use 2
    # transfinite points, leading to a width of 1 cell
    gmsh.model.geo.mesh.set_transfinite_curve(start_line, perp_nodes)

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
            line[i], normal_x, normal_y, cc_perpendicular_size
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
            line_length, cc_parallel_size
        )

        end_line = create_transfinite_cc_box(
            start_1, start_2, start_line,
            end_1, end_2,
            cc_loops, cc_line_surfaces, parallel_line_points, cc_lines,
            perp_nodes
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
    end_point_2, end_point_1 = get_extruded_points(line[-1], normal_x, normal_y, cc_perpendicular_size)

    # Create actual Gmsh points from the extrusion points
    end_1 = gmsh.model.geo.add_point(end_point_1[0], end_point_1[1], 0)
    end_2 = gmsh.model.geo.add_point(end_point_2[0], end_point_2[1], 0)

    # Compute the transfinite points of the parallel lines
    # We use the line from start_2 -> end_1 as our measuring stick, as
    # the difference between the parallel lines is likely very small
    line_length = distance(line[-1], line[-2])
    parallel_line_points = calculate_number_of_points(
        line_length, cc_parallel_size
    )

    # Create a transfinite surface out of the created CC box
    create_transfinite_cc_box(
        start_1, start_2, start_line,
        end_1, end_2,
        cc_loops, cc_line_surfaces, parallel_line_points,
        cc_lines, perp_nodes
    )