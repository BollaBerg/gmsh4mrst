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
