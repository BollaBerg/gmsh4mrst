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
    if 0 <= t <= 1 and 0 <= u <= 1:
        return [line_1_start[0] + t * delta_1[0], line_1_start[1] + t * delta_1[1]]
    else:
        return None


def get_perpendicular(delta_x, delta_y) -> 'tuple[float, float]':
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
