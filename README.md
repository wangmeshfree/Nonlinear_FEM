# MAE 5036: Advanced Computational Solid Mechanics II

This repository contains finite element method examples for MAE 5036. The runnable Python programs are in `src`; the Jupyter notebooks are step-by-step materials for classroom instruction.

## FEM Code

[src/FEM1D_linear_bar](src/FEM1D_linear_bar) | Runnable 1D linear-bar FEM program. See its [local guide](src/FEM1D_linear_bar/README.md). 

[src/FEM2D](src/FEM2D) | Shared runnable Q4 FEM program for the cantilever and plate-with-hole examples. See its [local guide](src/FEM2D/README.md).

Students: Please run the main.py function in this folder and try to understand how FEM solves 1 dimensional problem.


## Requirements

- Python 3
- NumPy
- Matplotlib
- ipykernel

Students are encouraged to install Python and learn to use it independently. You may use a different programming language for your work, provided you can implement and explain your solution. It is highly recommended to use VScode to write, run and debug the code for the future FEM practise. 

Step 1: Create a venv, the local environment for you to run the python code. 

Step 2: Install all the packages to this local environment
```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install ipykernel numpy matplotlib
```



The program creates a mesh-specific results folder under `src/FEM1D_linear_bar` containing nodal displacements plus displacement and stress plots.

## Classroom Notebooks

Open the notebooks in VS Code or Jupyter for the guided derivations, intermediate FEM calculations, plots, and convergence study used during class instruction.

[FEM1D_linear_bar.ipynb](notebooks/FEM1D_linear_bar.ipynb) | This is a classroom notebook for the 1D bar example. 

[FEM2D_linear_plate_with_hole.ipynb](notebooks/FEM2D_linear_plate_with_hole.ipynb) | This is a classroom notebook for the 2D plate with hole problem demonstration. 



## Author

Jiarui Wang, Zihan Wang

Copyright (c) 2026.
