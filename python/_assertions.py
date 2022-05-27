"""
    Helper files for internal use - Basic assertions.
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

def assert_column_in_dict(face_constraints: dict, column: str):
    if column not in face_constraints.keys():
        raise ValueError(
            f"Dictionary {dict} must contain column `{column}`!"
        )

def assert_columns_have_same_length(face_constraints: dict, col1: str, col2: str):
    if len(face_constraints.get(col1)) != len(face_constraints.get(col2)):
        raise ValueError(
            f"Column `{col1}` and `{col2}` must be the same length! "
            + f"len({col1}): {len(face_constraints.get(col1))}, len({col2}): {len(face_constraints.get(col2))}"
        )
