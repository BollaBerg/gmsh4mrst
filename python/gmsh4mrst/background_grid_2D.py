"""
    Create a 2D background mesh, without any embedding.
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
from typing import Any, Union, Iterable

import gmsh

from gmsh4mrst._arguments import (
    format_constraints, format_shape,
    format_meshing_algorithm, format_recombination_algorithm
)
from gmsh4mrst._geometry import split_at_intersections
from gmsh4mrst._gmsh import (
    create_threshold_field, create_circumference, create_fracture_point,
)

def background_grid_2D(
        cell_dimensions: float,
        *,
        shape: list = None,
        face_constraints: Union[
                    'list[list[Iterable]]',
                    'dict[str, float]',
                    'dict[str, Iterable]',
                    'dict[str, dict[str, float]]',
                    'dict[Any, dict[str, Iterable]]'] = None,
        face_constraint_factor: float = 1/3,
        min_FC_threshold_distance: float = 0.05,
        max_FC_threshold_distance: float = 0.2,
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
        min_CC_threshold_distance: float = 0.05,
        max_CC_threshold_distance: float = 0.2,
        CC_mesh_sampling: int = 100,
        mesh_algorithm: str = "Delaunay",
        recombination_algorithm: str = None,
        savename: str = None,
        run_frontend: bool = False
    ):
    """Create a 2D mesh, with user-supplied face- and cell constraints.

    Unlike delaunay_grid_2D, no constraints are actually embedded in the final
    result. Instead, only grid refinement is done. This is done due to the
    conversion from Delaunay triangulations to PEBI grids losing information
    about cell faces, leaving constraints messed up. By only refining the grid,
    we can instead handle constraints in MATLAB.

    This project was done as part of my Bachelor thesis during spring 2022, in
    order to help provide a secondary backend to the SINTEF-developed MATLAB
    module MRST, and its submodule UPR. For more information about MRST, see
    https://www.sintef.no/projectweb/mrst/. For more information about UPR, see
    https://www.sintef.no/projectweb/mrst/modules/upr/

    Args:
        cell_dimensions (float): Base dimensions of each cell.

        shape (Iterable[float] | Iterable[Iterable[float]] | dict[str, float],
                optional):
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
            If None, shape will be set to the square between (0, 0) and (1, 1).
            Defaults to None.
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

        min_FC_threshold_distance (float, optional): Distance from face
            constraints where cell dimensions will start increasing. Defaults
            to 0.05.

        max_FC_threshold_distance (float, optional): Distance from face
            constraints where cell dimensions will be back to their default
            (max) value, i.e. the supplied argument cell_dimensions. Defaults
            to 0.2.

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
            around the cell constraint lines, as compared to the supplied
            cell_dimensions. Overrides cell_constraint_factor for lines. If set
            to None, cell_constraint_factor will be used for lines. Defaults
            to None.
        
        cell_constraint_point_factor (float, optional): The size used for cells
            around cell constraint points, as compared to the supplied
            cell_dimensions. Overrides cell_constraint_factor for points. If
            set to None, cell_constraint_factor will be used for points.
            Defaults to None.

        min_CC_threshold_distance (float, optional): Distance from cell
            constraints where cell dimensions will start increasing. Defaults
            to 0.05.

        max_CC_threshold_distance (float, optional): Distance from cell
            constraints where cell dimensions will be back to their default
            (max) value, i.e. the supplied argument cell_dimensions. Defaults
            to 0.2.

        CC_mesh_sampling (int, optional): The number of points along the
            cell constraints should be sampled to calculate the threshold
            distances. Defaults to 100.

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
            be saved. Defaults to None.

        run_frontend (bool, optional): Set to True in order to run the Gmsh
            frontend and show the created mesh. Defaults to False.
    """
    # Handle default values
    if shape is None:
        shape = [(0, 0), (1, 0), (1, 1), (0, 1)]
    if face_constraints is None:
        face_constraints = []
    if min_intersection_distance is None:
        min_intersection_distance = min_FC_threshold_distance
    if max_intersection_distance is None:
        max_intersection_distance = max_FC_threshold_distance
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

    # Split all constraints at their intersection
    face_constraints, cell_constraints, intersection_IDs = split_at_intersections(
        face_constraints, cell_constraints
    )
    
    gmsh.initialize()
    gmsh.model.add("gmsh4mrst")

    # Create corners
    corners = [
        gmsh.model.geo.add_point(point[0], point[1], 0) for point in shape
    ]

    # Create circumference
    circumference = create_circumference(corners)

    ##### CREATE FRACTURES (FACE CONSTRAINTS) ######
    fracture_points = []
    fractures = []
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
        fracture_start = create_fracture_point(line[0], intersection_IDs)
        for i in range(1, len(line)):
            fracture_end = create_fracture_point(line[i], intersection_IDs)
            fractures.append(
                gmsh.model.geo.add_line(fracture_start, fracture_end)
            )
            fracture_start = fracture_end

    # Create cell constraints\
    cc_lines = []               # Holds lines - used when refining mesh
    cc_points = []              # Holds single points - used when refining mesh
    cc_point_size = cell_constraint_point_factor * cell_dimensions  # Cell size around points
    cc_line_size = cell_constraint_line_factor * cell_dimensions    # Cell size around lines
    for line in cell_constraints:
        if len(line) == 1:
            # line is a single point
            cc_points.append(
                gmsh.model.geo.add_point(line[0][0], line[0][1], 0)
            )
        
        else:
            # line has at least 1 line segment
            line_start = create_fracture_point(line[0], intersection_IDs)
            for i in range(1, len(line)):
                line_end = create_fracture_point(line[i], intersection_IDs)
                cc_lines.append(
                    gmsh.model.geo.add_line(line_start, line_end)
                )
                line_start = line_end

    # Create curve loop of circumference
    circumference_loop = gmsh.model.geo.add_curve_loop(circumference)
    
    # Define surface from circumference
    surface = gmsh.model.geo.add_plane_surface([circumference_loop])

    # Synchronize to prepare for embedding fracture face_constraints
    gmsh.model.geo.synchronize()

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
        min_distance=min_FC_threshold_distance,
        max_distance=max_FC_threshold_distance
    )

    # Add field for intersection points
    intersection_points = list(intersection_IDs.values())
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

    # Refine mesh around cell constraint lines
    cell_line_threshold = create_threshold_field(
        None,
        curve_list=cc_lines,
        sampling=CC_mesh_sampling,
        min_size=cc_line_size,
        max_size=cell_dimensions,
        min_distance=min_CC_threshold_distance,
        max_distance=max_CC_threshold_distance,
    )

    # Refine mesh around cell constraint points
    cell_point_threshold = create_threshold_field(
        cc_points,
        curve_list=None,
        sampling=CC_mesh_sampling,
        min_size=cc_point_size,
        max_size=cell_dimensions,
        min_distance=min_CC_threshold_distance,
        max_distance=max_CC_threshold_distance,
    )

    # We use the minimum of all fields as our background mesh
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [
        fracture_threshold, intersection_threshold,
        cell_line_threshold, cell_point_threshold
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
    background_grid_2D(
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
        mesh_algorithm="DelQuad",
        savename=None,
        run_frontend=True)