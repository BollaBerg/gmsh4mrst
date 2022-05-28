"""
    Create a base mesh, optimized for converting to PEBI later.
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
    create_cell_constraint_point, create_cell_constraint_line
)

def pebi_base_2D(
        cell_dimensions: float,
        *,
        shape: list = None,
        face_constraints: Union[
                    'list[list[Iterable]]',
                    'dict[str, float]',
                    'dict[str, Iterable]',
                    'dict[str, dict[str, float]]',
                    'dict[Any, dict[str, Iterable]]'] = None,
        face_constraint_factor: float = 1/4,
        face_constraint_parallel_factor: float = None,
        face_constraint_perpendicular_factor: float = None,
        face_constraint_perpendicular_cells: int = 3,
        face_constraint_point_factor: float = None,
        face_constraint_refinement_factor: float = None,
        min_FC_threshold_distance: float = 0.05,
        max_FC_threshold_distance: float = 0.2,
        face_mesh_sampling: int = 100,
        cell_constraints: Union[
                    'list[list[Iterable]]',
                    'dict[str, float]',
                    'dict[str, Iterable]',
                    'dict[str, dict[str, float]]',
                    'dict[Any, dict[str, Iterable]]'] = None,
        cell_constraint_factor: float = 1/4,
        cell_constraint_parallel_factor: float = None,
        cell_constraint_perpendicular_factor: float = None,
        cell_constraint_perpendicular_cells: int = 2,
        cell_constraint_point_factor: float = None,
        cell_constraint_refinement_factor: float = None,
        min_CC_threshold_distance: float = 0.05,
        max_CC_threshold_distance: float = 0.2,
        CC_mesh_sampling: int = 100,
        mesh_algorithm: str = "Delaunay",
        recombination_algorithm: str = None,
        savename: str = None,
        run_frontend: bool = False
    ):
    """Create a 2D base for converting to a PEBI grid, with user-supplied face-
    and cell constraints.

    Face constraints are traced a surrounding transfinite mesh, with three
    nodes on each short side. This ensures that the PEBI faces will align with
    the constraint as well as possible.
    Cell constraints are traced by a surrounding transfinite mesh with four
    nodes on each short site. This ensures that the PEBI cells will align with
    the constraint as well as possible. 

    NOTE: This is meant as a base for later conversion to PEBI grid, primarily
    in UPR. As such, face constraints and cell constraints are swapped and
    reinforced, as an attempt at making the conversion as accurate as possible.

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
            Constraints for the PEBI faces. Each constraint is the coordinate
            of a surface-trace. The surface is assumed to be linear between 
            the coordinates. The resulting mesh will place sites such that 
            the surfaces are traced by faces of the grid when converted to PEBI.
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

        face_constraint_parallel_factor (float, optional): The size used
            for the length of the cells within the face constraints, as
            compared to the supplied cell_dimensions. Overrides
            face_constraint_factor for lines. If set to None,
            face_constraint_factor will be used for length. Defaults to None.
        
        face_constraint_perpendicular_factor (float, optional): The size used
            for the width of the cells within the face constraints, as
            compared to the supplied face_dimensions. Overrides
            face_constraint_factor for lines. If set to None,
            face_constraint_factor will be used for width. Defaults to None.
        
        face_constraint_perpendicular_cells (int, optional): The number of 
            transfinite cells to use across the face constraints. Should be
            odd if the Delaunay grid is used directly, and even if the grid is
            converted to a PEBI grid before use. More cells may make a "buffer
            zone" around the constraints. Defaults to 3.
            
        face_constraint_point_factor (float, optional): The size used for cells
            around face constraint points, as compared to the supplied
            cell_dimensions. Overrides face_constraint_factor for points. If
            set to None, face_constraint_factor will be used for points.
            Defaults to None.

        face_constraint_refinement_factor (float, optional): The cell size in
            the refinement along the face constraint lines. If set to None,
            cell_dimensions will be used and no refinement will be done.
            Defaults to None.

        min_FC_threshold_distance (float, optional): Distance from face
            constraints where cell dimensions will start increasing. Defaults
            to 0.05.

        max_FC_threshold_distance (float, optional): Distance from face
            constraints where cell dimensions will be back to their default
            (max) value, i.e. the supplied argument cell_dimensions. Defaults
            to 0.2.

        fracture_mesh_sampling (int, optional): The number of points along the
            face constraints should be sampled to calculate the threshold
            distances. Defaults to 100.

        cell_constraints (list[Iterable] | dict[str, float] | dict[str, Iterable]
                | dict[str, dict[str, float]] | dict[str, dict[str, Iterable]],
                optional):
            Constraints for the cell centroids. Each constraint is a line,
            which the method will attempt to place in the center of the returned
            grid when converted to PEBI. The lines are assumed to be linear
            between the coordinates. If a constraint is only one coordinate,
            the line is treated as a point constraint.
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

        cell_constraint_parallel_factor (float, optional): The size used
            for the length of the cells within the cell constraints, as
            compared to the supplied cell_dimensions. Overrides
            cell_constraint_factor for lines. If set to None,
            cell_constraint_factor will be used for length. Defaults to None.
        
        cell_constraint_perpendicular_factor (float, optional): The size used
            for the width of the cells within the cell constraints, as
            compared to the supplied cell_dimensions. Overrides
            cell_constraint_factor for lines. If set to None,
            cell_constraint_factor will be used for width. Defaults to None.
        
        cell_constraint_perpendicular_cells (int, optional): The number of 
            transfinite cells to use across the cell constraints. Should be
            even if the Delaunay grid is used directly, and odd if the grid is
            converted to a PEBI grid before use. More cells may make a "buffer
            zone" around the constraints. Defaults to 2.
            
        cell_constraint_point_factor (float, optional): The size used for cells
            around cell constraint points, as compared to the supplied
            cell_dimensions. Overrides cell_constraint_factor for points. If
            set to None, cell_constraint_factor will be used for points.
            Defaults to None.

        cell_constraint_refinement_factor (float, optional): The cell size in
            the refinement along the cell constraint lines. If set to None,
            cell_dimensions will be used and no refinement will be done.
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
            be saved. The MATLAB functions assume that the file will be saved
            as "TEMP_Gmsh_MRST.m". Defaults to "TEMP_Gmsh_MRST.m".

        run_frontend (bool, optional): Set to True in order to run the Gmsh
            frontend and show the created mesh. Defaults to False.
    """
    # Handle default values
    if shape is None:
        shape = [(0, 0), (1, 0), (1, 1), (0, 1)]
    if face_constraints is None:
        face_constraints = []
    if face_constraint_parallel_factor is None:
        face_constraint_parallel_factor = face_constraint_factor
    if face_constraint_perpendicular_factor is None:
        face_constraint_perpendicular_factor = face_constraint_factor
    if face_constraint_point_factor is None:
        face_constraint_point_factor = face_constraint_factor
    if face_constraint_refinement_factor is None:
        face_constraint_refinement_factor = 1
    if cell_constraints is None:
        cell_constraints = []
    if cell_constraint_parallel_factor is None:
        cell_constraint_parallel_factor = cell_constraint_factor
    if cell_constraint_perpendicular_factor is None:
        cell_constraint_perpendicular_factor = cell_constraint_factor
    if cell_constraint_point_factor is None:
        cell_constraint_point_factor = cell_constraint_factor
    if cell_constraint_refinement_factor is None:
        cell_constraint_refinement_factor = 1

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
    face_loops = []               # Holds line loops - used when making surface
    face_lines = []               # Holds lines - used when refining mesh
    face_point_surfaces = []      # Holds surfaces around CC points
    face_line_surfaces = []       # Holds surfaces around CC lines
    face_point_size = face_constraint_point_factor * cell_dimensions                    # Cell size around points
    face_parallel_size = face_constraint_parallel_factor * cell_dimensions              # Cell size along lines
    face_perpendicular_size = face_constraint_perpendicular_factor * cell_dimensions    # Cell size across lines
    for line in face_constraints:
        if len(line) == 1:
            # line is a single point
            create_cell_constraint_point(
                line[0], face_point_size, face_loops, face_point_surfaces, face_lines,
                nodes=face_constraint_perpendicular_cells
            )
        
        else:
            # line has at least 1 line segment
            create_cell_constraint_line(
                line, face_parallel_size, face_perpendicular_size, face_loops,
                face_line_surfaces, face_lines, perp_nodes=face_constraint_perpendicular_cells
            )

    # Create cell constraints
    cc_loops = []               # Holds line loops - used when making surface
    cc_lines = []               # Holds lines - used when refining mesh
    cc_point_surfaces = []      # Holds surfaces around CC points
    cc_line_surfaces = []       # Holds surfaces around CC lines
    cc_point_size = cell_constraint_point_factor * cell_dimensions  # Cell size around points
    cc_parallel_size = cell_constraint_parallel_factor * cell_dimensions    # Cell size along lines
    cc_perpendicular_size = cell_constraint_perpendicular_factor * cell_dimensions    # Cell size across lines
    for line in cell_constraints:
        if len(line) == 1:
            # line is a single point
            create_cell_constraint_point(
                line[0], cc_point_size, cc_loops, cc_point_surfaces, cc_lines,
                nodes=cell_constraint_perpendicular_cells
            )
        
        else:
            # line has at least 1 line segment
            create_cell_constraint_line(
                line, cc_parallel_size, cc_perpendicular_size, cc_loops,
                cc_line_surfaces, cc_lines, perp_nodes=cell_constraint_perpendicular_cells
            )

    # Create curve loop of circumference
    circumference_loop = gmsh.model.geo.add_curve_loop(circumference)
    
    # Define surface from circumference
    # We remove the loops created from the constraints
    surface = gmsh.model.geo.add_plane_surface([circumference_loop, *face_loops, *cc_loops])

    # Synchronize to prepare for embedding fracture face_constraints
    gmsh.model.geo.synchronize()

    # Refine mesh around fractures. We want a finer mesh around the fractures,
    # to ensure as close fit of faces as possible. In addition, this helps
    # handle intersecting fractures
    # This is gotten from
    # https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t10.py

    face_threshold = create_threshold_field(
        point_list=None,
        curve_list=face_lines,
        sampling=face_mesh_sampling,
        min_size=face_constraint_refinement_factor * cell_dimensions,
        max_size=cell_dimensions,
        min_distance=min_FC_threshold_distance,
        max_distance=max_FC_threshold_distance
    )


    # Refine mesh around cell constraints as well
    cell_threshold = create_threshold_field(
        None,
        curve_list=cc_lines,
        sampling=CC_mesh_sampling,
        min_size=cell_constraint_refinement_factor * cell_dimensions,
        max_size=cell_dimensions,
        min_distance=min_CC_threshold_distance,
        max_distance=max_CC_threshold_distance,
    )

    # We use the minimum of all fields as our background mesh
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [
        face_threshold, cell_threshold
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
    pebi_base_2D(
        0.2, 
        face_constraints=[
            [(0.25, 0.25), (0.4, 0.5), (0.7, 0.7)],
            [(0.8, 0.1), (0.9, 0.2)]
        ],
        cell_constraints=[
            [(0.1, 0.1)],
            [(0.5, 0.8), (0.6, 0.7), (0.7, 0.8), (0.9, 0.6)],
        ],
        shape=[
            (0, 0), (0.5, 0.2), (1, 0), (1, 1), (0, 1)
        ],
        savename=None,
        run_frontend=True)