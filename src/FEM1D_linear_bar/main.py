# ******************************************************************************
# *                                                                            *
# *        MAE 5036 - Advanced Computational Solid Mechanics II                *
# *                                                                            *
# *                  1D Bar Finite Element Analysis                            *
# *                             JIARUI WANG                                    *
# *                         Copyright (c) 2026                                 *
# *                                                                            *
# ******************************************************************************

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import fem1d_lib as fem1d

E = 100.0
A = 1.0
L = 10.0
P = 5.0


# 1. Choose the mesh density and integration rule.
number_element_in_one_direction = 2
number_gauss_points_in_one_direction = 2

output_directory = (
    Path(__file__).resolve().parent
    / f"Results_1d_bar_{number_element_in_one_direction}_elements"
)
output_directory.mkdir(exist_ok=True)



# 2. Generate nodes, connectivity, and boundary-condition data.
total_nodes, coords, elements, fixed_nodes, fixed_values, nodal_forces = (
    fem1d.get_mesh_bar(
        number_element_in_one_direction,
        L,
        P,
    )
)
total_element = len(elements)

# 3. Generate Gauss points and weights on [-1, 1].
gauss_points, gauss_weights = fem1d.get_gauss_integration_points(
    number_gauss_points_in_one_direction,
)

# 4. Assemble stiffness and the equivalent nodal distributed loads.
K, F = fem1d.get_global_stiffness_and_force(
    total_nodes,
    coords,
    total_element,
    elements,
    E,
    A,
    gauss_points,
    gauss_weights,
    fem1d.get_body_force,
)
# 5. Add the prescribed concentrated force at the right end.
F = fem1d.apply_traction_boundary_conditions(F, nodal_forces)

# 6. Apply prescribed displacements (this modifies K and F).
K, F = fem1d.apply_essential_boundary_conditions(
    K,
    F,
    fixed_nodes,
    fixed_values,
)

# 7. Solve for the nodal displacements.
U = np.linalg.solve(K, F)

print("Number of nodes:", total_nodes)
print("Number of elements:", total_element)
print("Node coordinates:", coords)
print("FEM nodal displacements:", U)

displacement_data = np.column_stack((coords, U))
np.savetxt(
    output_directory / "nodal_displacements.txt",
    displacement_data,
    header="x  FEM_displacement",
    fmt="%.8e",
)


# 8. Visualization 
x_exact = np.linspace(0.0, L, 300)

sample_coordinates, sample_displacements, _, sample_stresses = (
    fem1d.compute_element_displacement_strain_and_stress(
        coords,
        elements,
        U,
        E,
    )
)

fem1d.plot_element_variables(
    sample_coordinates,
    sample_displacements,
    x_exact,
    fem1d.exact_displacement(x_exact, E, A, L, P),
    coords,
    "FEM displacement",
    "Elementwise Axial Displacement",
    output_directory / "displacement.png",
)

fem1d.plot_element_variables(
    sample_coordinates,
    sample_stresses,
    x_exact,
    fem1d.exact_stress(x_exact, L, P, A),
    coords,
    "FEM stress",
    "Elementwise Axial Stress (No Nodal Averaging)",
    output_directory / "stress.png",
)

