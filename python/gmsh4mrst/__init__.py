"""
    Automatically create Gmsh meshes

    These can be used directly in a Python program, but are mainly designed to
    work with the MATLAB module MRST.  
    
    FUNCTIONS:
        background_grid_2D
        delaunay_grid_2D
        pebi_base_2D

    LICENSE:
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

from .background_grid_2D import background_grid_2D
from .delaunay_grid_2D import delaunay_grid_2D
from .pebi_base_2D import pebi_base_2D