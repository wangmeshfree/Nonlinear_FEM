# 2D FEM Examples

This folder contains one shared Q4 finite-element implementation for both the plate-with-hole and cantilever examples. The solver functions are in `fem2d_lib.py`; `main.py` selects the problem, input mesh, parameters, and named boundary sets.

## Files

| File | Purpose |
| --- | --- |
| `main.py` | Shared executable entry point. Change `active_problem` to select the example. |
| `fem2d_lib.py` | Mesh reader, Q4 shape functions, integration, assembly, boundary conditions, postprocessing, and error calculation. |
| `Inputs/` | LS-DYNA keyword meshes and Gmsh mesh files. |

## Select a problem

In `main.py`:

```python
active_problem = 'Cantilever_beam'
```

or:

```python
active_problem = 'Plate_with_hole'
```

Each entry in `problem_settings` supplies the mesh filename, analytical-solution parameters, and the physical-group names used for Dirichlet and traction boundaries.

The current constitutive matrix is the plane-strain matrix used in the classroom notebook. The two examples share the same Q4 assembly code; their differences are the analytical field, mesh, parameters, and boundary sets.

## Run

From this directory:

```sh
python main.py
```

The environment must contain NumPy, Matplotlib, and ipykernel. Gmsh is needed when regenerating a mesh from a `.geo` file.
