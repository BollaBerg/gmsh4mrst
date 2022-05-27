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
