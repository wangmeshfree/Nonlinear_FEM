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

The plate-with-hole example uses plane strain. The cantilever example uses plane stress and applies the analytical two-dimensional displacement field on the left edge using the `(99, 99)` Dirichlet flags. In the boundary-condition flags, `0` means free, `1` means zero displacement, and `99` means analytical displacement. The two examples share the same Q4 assembly code; their differences are the analytical field, mesh, parameters, constitutive matrix, and named boundary sets.

`main.py` contains the same step-by-step workflow as the 1D example: choose a configuration, create the constitutive matrix, read the mesh, assemble and solve one mesh, save the solution, then run the configured mesh-convergence study. `fem2d_lib.py` contains reusable mesh, Q4, assembly, boundary-condition, plotting, and error functions.

## Run

From this directory:

```sh
python main.py
```

The environment must contain NumPy, Matplotlib, and ipykernel. Gmsh is needed when regenerating a mesh from a `.geo` file.

## Results

Running `main.py` creates a folder beside it, for example:

```text
Results_2d_Cantilever_beam_model_cantilever_level_2/
├── nodal_displacements.txt
├── nodal_stresses.txt
├── mesh.png
├── displacement.png
├── stress.png
├── convergence.txt
└── convergence.png
```

The text files contain the numerical values. The PNG files contain the mesh, displacement, stress, and convergence plots. The folder name changes automatically with the selected problem and main mesh.
