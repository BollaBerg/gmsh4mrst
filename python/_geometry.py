"""
    Helper files for internal use - Geometry related.
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

from copy import deepcopy
from math import atan2, sqrt, ceil


def find_intersection(line_1_start, line_1_end, line_2_start, line_2_end):
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
    # Set to <= to also detect corner points. Set to < to only detect internals
    if 0 < t < 1 and 0 < u < 1:
        return (line_1_start[0] + t * delta_1[0], line_1_start[1] + t * delta_1[1])
    else:
        return None


def get_perpendicular(delta_x, delta_y) -> 'tuple[float, float]':
    """Return the normal vector of delta_x, delta_y pointing 'left'"""
    return -delta_y, delta_x


def get_extruded_points(
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

def line_bends_towards_right(start_point, mid_point, end_point) -> bool:
    angle_start_mid = atan2(
        mid_point[1] - start_point[1], mid_point[0] - start_point[0]
    )
    angle_start_end = atan2(
        end_point[1] - start_point[1], end_point[0] - start_point[0]
    )
    return angle_start_mid > angle_start_end


def distance(point_1, point_2) -> float:
    """Compute the distance between two points, using Euclidean distance"""
    return sqrt((point_2[1] - point_1[1])**2 + (point_2[0] - point_1[0])**2)


def calculate_number_of_points(line_length, wanted_segment_size):
    return ceil(line_length / wanted_segment_size) + 1


def get_midpoint(point_1: 'tuple[float, float]',
                 point_2: 'tuple[float, float]') -> 'tuple[float, float]':
    return (
        0.5 * (point_1[0] + point_2[0]),
        0.5 * (point_1[1] + point_2[1])
    )

def split_at_intersections(
        face_constraints: list, cell_constraints: list
        ) -> 'tuple[list, list, dict]':
    """Split a formatted list of constraints at their intersections

    Args:
        constraints (list): Formatted list of contraints. Must have shape:
            [
                [(x11, y11), (x12, y12), ...],
                [(x21, y21), (x22, y22), ...],
                ...
            ]

    Returns:
        list: Formatted face constraints, but with every constraint broken
            into non-intersecting (except possibly at the end) line segments
        list: Formatted cell constraints, but with every constraint broken
            into non-intersecting (except possibly at the end) line segments
        dict: Dictionary of intersection IDs, used to avoid creating multiple
            Gmsh points in each intersection
    """
    # We save intersections in a dict
    # That way, we avoid creating multiple Gmsh point in each intersection
    intersection_IDs = dict()

    # Avoid changing the inputs
    face_const = deepcopy(face_constraints)
    cell_const = deepcopy(cell_constraints)

    def _check_line(line, other_constraints):
        if len(line) < 2:
            # line is a point
            return
        segment_ID = 0
        while segment_ID < len(line) - 1:
            start = line[segment_ID]         # start = (x11, y11)
            end = line[segment_ID+1]
            # Check for intersections in own line
            # We only need to check "forwards", as we have already looked at
            # the combinations behind us
            # We naturally don't check intersections with the current segment
            other_ID = segment_ID + 1
            while other_ID < len(line) - 1:
                other_start = line[other_ID]
                other_end = line[other_ID + 1]

                intersection = find_intersection(start, end, other_start, other_end)
                if intersection is not None:
                    intersection_IDs[intersection] = None
                    # Insert intersection into both segments of the line
                    line.insert(segment_ID + 1, intersection)
                    line.insert(other_ID + 2, intersection)

                    # Get new segment end 
                    end = line[segment_ID + 1]
                other_ID += 1

            # check for intersections in the other constraints
            for other_line in other_constraints:
                other_ID = 0
                while other_ID < len(other_line) - 1:
                    other_start = other_line[other_ID]
                    other_end = other_line[other_ID + 1]

                    intersection = find_intersection(start, end, other_start, other_end)
                    if intersection is not None:
                        intersection_IDs[intersection] = None
                        # Insert intersection into both segments of the line
                        line.insert(segment_ID + 1, intersection)
                        other_line.insert(other_ID + 1, intersection)

                        # Get new segment end 
                        end = line[segment_ID + 1]
                    other_ID += 1
            segment_ID += 1

    # Go through face constraints first
    for line_id, line in enumerate(face_const):
        _check_line(line, face_const[line_id + 1: ] + cell_const)

    # Now check cell constraints against each other
    for line_id, line in enumerate(cell_const):
        _check_line(line, cell_const[line_id + 1: ])
    return face_const, cell_const, intersection_IDs
