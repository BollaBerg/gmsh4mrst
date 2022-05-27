"""
    Helper files for internal use - Argument formatting.
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

from array import array
from typing import Iterable, Union

from gmsh4mrst._assertions import (
    assert_column_in_dict, assert_columns_have_same_length
)

def format_constraints(constraints) -> 'list[list[tuple[float, float]]]':
    """Check that constraints has one of the legal forms, and format it

    Args:
        constraints (list[tuple]
                     list[list[tuple]]
                    | dict[str, float]
                    | dict[str, Iterable]
                    | dict[str, dict[str, float]]
                    | dict[str, dict[str, Iterable]]
            ): constraints that should be formatted

    Returns:
        list: The constraints in the form
            [
                [(x11, y11), (x12, y12), ...],
                [(x21, y21), (x22, y22), ...],
                ...
            ]
    """
    if isinstance(constraints, dict):
        if len(constraints) == 0:
            # constraints is an empty list -> No fractures to handle,
            # simply create a single mesh within size
            output = []
        elif isinstance(list(constraints.values())[0], array):
            # If data is sent from MATLAB, then constraints is a dict
            # If there is only one line, then we accept it as a 1D dict, i.e.
            # one can write
            #   constraints.x = [x1 ...];
            #   constraints.y = [y1 ...];
            # Then constraints is a dict in the shape {x: array, y: array}
            output = [_check_array_dict_and_return_line_list(constraints)]

        elif isinstance(list(constraints.values())[0], dict):
            # If data is sent from MATLAB, then constraints is a dict of dicts
            # If there are more than one line, then it must be a 2D dict, i.e.
            # one can write
            #   constraints.line1.x = [x11 ...];
            #   constraints.line1.y = [y11 ...];
            #   constraints.line2.x = [x21 ...];
            #   constraints.line2.y = [y21 ...];
            # Each element in constraints.values is a dict of {x: array, y: array}
            output = []
            for constraint in constraints.values():
                output.append(
                    _check_array_dict_and_return_line_list(constraint)
                )
    elif isinstance(constraints, Iterable):
        if len(constraints) == 0:
            # constraints is an empty list -> No fractures to handle,
            # simply create a single mesh within size
            output = []
        elif len(constraints) == 2 and isinstance(constraints[0], float):
            # constraints has form [x, y], likely one point
            output = [(constraints[0], constraints[1])]
        elif isinstance(constraints[0], tuple):
            # constraints is an iterable of tuples or a single point
            # -> constraints is a single list of points
            # -> check that every point contains x- and y-value
            output = []
            for constraint in constraints:
                if not isinstance(constraint, Iterable) or len(constraint) < 2:
                    raise ValueError(
                        f"Constraint {constraint} must be in the form (x, y)"
                    )
                output.append((constraint[0], constraint[1]))
        elif isinstance(constraints[0], Iterable):
            # constraints is a list lines
            # -> check that each point has x- and y-value
            output = []
            for line in constraints:
                output_line = []
                for point in line:
                    if not isinstance(point, Iterable) or len(point) < 2:
                        raise ValueError(
                            f"Constraint {point} must be in the form (x, y)"
                        )
                    output_line.append((point[0], point[1]))
                output.append(output_line)
    else:
        raise ValueError(
            f"Received unknown constraint {constraints}"
        )
    return output


def format_shape(shape) -> 'list[tuple[float, float]]':
    """Check that shape has one of the legal forms, and format it

    Args:
        shape (Iterable[float]
               | Iterable[Iterable[float]]
               | dict[str, float]
            ): Shape to format

    Raises:
        ValueError: Raised if shape has illegal type or shape

    Returns:
        The shape in the form [(x1, y1), (x2, y2), ...]
    """
    if isinstance(shape, dict):
        # Shape is a dict, likely from MATLAB
        assert_column_in_dict(shape, "x")
        assert_column_in_dict(shape, "y")
        assert_columns_have_same_length(shape, "x", "y")
        output = [
            (x, y) for x, y in zip(shape.get("x"), shape.get("y"))
        ]
    elif isinstance(shape, Iterable):
        if len(shape) < 2:
            raise ValueError(
                f"Shape must have length >= 2. Current length: {len(shape)}"
            )
        elif len(shape) == 2:
            # Shape is the size of the domain, starting at (0, 0)
            output = [
                (0, 0), (0, shape[1]), (shape[0], shape[1]), (shape[0], 0)
            ]
        else:
            # output is a list of points making up the domain
            # We check that each point consists of two actual points
            output = []
            for point in shape:
                if not isinstance(point, Iterable) or len(point) < 2:
                    raise ValueError(
                        f"Point {point} must be in the form (x, y)"
                    )
                output.append((point[0], point[1]))
    else:
        raise ValueError(
            f"Shape {shape} (type {type(shape)}) is not supported for argument `shape`!"
        )
    return output


def format_meshing_algorithm(mesh_algorithm: Union[str, int]) -> int:
    """Check that mesh_algorithm has legal value, and format it
    
    According to 
    https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t10.py
    Frontal-Delaunay (6) usually leads to the highest quality meshes, but
    Delaunay (5) handles complex mesh sizes better - especially size fields
    with large element size gradients. We therefore default to Delaunay.
    For quad-shaped grids, "Frontal-Delaunay for Quads" (8) may be beneficial.

    Args:
        mesh_algorithm (str | int): Input to format

    Raises:
        ValueError: Raise if mesh_algorithm has illegal type or value

    Returns:
        int: Gmsh-ID of the meshing algorithm
    """
    mapping = {
        "meshadapt": 1,
        "automatic": 2,
        "delaunay": 5,
        "frontal": 6,
        "bamg": 7,
        "delquad": 8
    }
    if mesh_algorithm is None:
        # Use Delaunay as default
        return 5
    elif isinstance(mesh_algorithm, int):
        # Check if mesh_algorithm is a legal value. If so, use it
        if mesh_algorithm not in mapping.values():
            raise ValueError(
                "mesh_algorithm must be a legal value. Current value: "
                + f"{mesh_algorithm}. Legal values: {mapping}"
            )
        return mesh_algorithm
    elif isinstance(mesh_algorithm, str):
        if mesh_algorithm.lower() not in mapping.keys():
            raise ValueError(
                "mesh_algorithm must be a legal value. Current value: "
                + f"{mesh_algorithm}. Legal values: {mapping}"
            )
        return mapping.get(mesh_algorithm.lower())
    else:
        raise ValueError(
            f"Unsupported type for mesh_algorithm: {type(mesh_algorithm)}. "
            + "mesh_algorithm must be an int or a string."
        )


def format_recombination_algorithm(algorithm: Union[str, int]) -> int:
    """Check that recombination_algorithm has legal value, and format it
    
    Args:
        algorithm (str | int): Input to format

    Raises:
        ValueError: Raise if algorithm has illegal type or value

    Returns:
        int: Gmsh-ID of the recombination algorithm
    """
    mapping = {
        "simple": 0,
        "blossom": 1,
        "simplefull": 2,
        "blossomfull": 3
    }
    if algorithm is None:
        # Default is None
        return None
    elif isinstance(algorithm, int):
        if algorithm not in mapping.values():
            raise ValueError(
                "recombination_algorithm must be a legal value. Current value: "
                + f"{algorithm}. Legal values: {mapping}"
            )
        return algorithm
    elif isinstance(algorithm, str):
        if algorithm.lower() not in mapping.keys():
            raise ValueError(
                "recombination_algorithm must be a legal value. Current value: "
                + f"{algorithm}. Legal values: {mapping}"
            )
        return mapping.get(algorithm.lower())
    else:
        raise ValueError(
            f"Unsupported type for recombination_algorithm: {type(algorithm)}. "
            + "recombination_algorithm must be an int or a string."
        )


def format_constraint_factor(factor: Union[float, Iterable],
                             constraint_length: int) -> list:
    if isinstance(factor, float):
        return [factor] * constraint_length
    elif isinstance(factor, Iterable) and len(factor) == constraint_length:
        return factor
    elif isinstance(factor, Iterable):
        raise ValueError(
            f"Factor {factor} has wrong length. Expected length: "
            + f"{constraint_length}. Received: {len(factor)}"
        )
    else:
        raise ValueError(
            f"Factor ({factor}) of unexpected type {type(factor)} received. "
            + f"Expected type float or Iterable with length: {constraint_length}"
        )
    


def _check_array_dict_and_return_line_list(face_constraints: dict):
    # Check that face_constraints contains column x and y
    assert_column_in_dict(face_constraints, "x")
    assert_column_in_dict(face_constraints, "y")
    
    # Handle actual content
    if isinstance(face_constraints.get("x"), float):
        # If the 'line' is a single point, then x and y are floats,
        # and must be handled explicitly
        return [(face_constraints.get("x"), face_constraints.get("y"))]
    else:
        # Assume 'line' is a list/tuple/Iterable of floats
        assert_columns_have_same_length(face_constraints, "x", "y")
        return list(zip(face_constraints.get("x"), face_constraints.get("y")))