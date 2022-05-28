# gmsh4mrst

Automatically create Gmsh meshes for use in MRST.

## Building and pushing PyPi package
All package building, packaging and pushing to PyPi follows the instructions on [packaging.python.org](https://packaging.python.org/en/latest/tutorials/packaging-projects/). Here is a brief overview of commands needed:

Start by ensuring the version number has been updated (follow [Semantic Versioning](https://semver.org/))! This is done in `setup.cfg`.

Ensure the necessary packages are installed
```bash
pip install --upgrade build     # Used for building
pip install --upgrade twine     # Used for uploading
```

Build the project
```bash
python3 -m build
```
There should now be two files in the `dist/` directory:
```bash
dist/
    gmsh4mrst-VERSION-py3-none-any.whl
    gmsh4mrst-VERSION.tar.gz
```

Upload the project to PyPi
```bash
python3 -m twine upload dist/*
```