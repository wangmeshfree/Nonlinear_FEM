"""2D Q4 FEM driver: select a problem, solve, plot, and save results."""

from pathlib import Path
import numpy as np
import fem2d_lib as fem2d
from fem2d_problem_lib import problem_settings


# Problem settings live in fem2d_problem_lib.py; see PROBLEM_SETTINGS_EXAMPLES.md for templates.
# 1. Choose the problem and material.

# active_problem = 'Beam_with_tip_load'

active_problem = 'Cantilever_beam'
settings = problem_settings[active_problem]


E, nu = settings['material_properties']

number_gauss_points= settings['integration_order']

plane_problem_type = settings['2d_problem_type']

has_analytical = 'solution' in settings

if has_analytical:
    exact = settings['solution'](
        *settings['model_characteristics'], E, nu, plane_problem_type
    )
    body_force = exact.body_force
else:
    exact = None
    body_force_value = np.array(settings.get('body_force', (0.0, 0.0)))
    def body_force(x, y):
        return body_force_value
    

input_directory = Path(__file__).resolve().parent / 'Inputs'
output_directory = (
    Path(__file__).resolve().parent
    / f"Results_2d_{active_problem}_{Path(settings['mesh_file']).stem}"
)
output_directory.mkdir(exist_ok=True)


# 2. Build the constitutive matrix and Gauss rule.
D_mat = fem2d.get_elasticity_matrix(E, nu, plane_problem_type)
gauss_points, gauss_weights = fem2d.get_gauss_integration_points(
    number_gauss_points
)


# 3. Read the selected LS-DYNA keyword mesh and its named boundary sets.
total_nodes, coords, total_element, elements, boundary_connect, nodeset_list = (
    fem2d.read_lsdyna_k(input_directory / settings['mesh_file'])
)
fix_nodes = fem2d.get_fixnodes_from_sets(
    nodeset_list, settings['dirichlet'], coords, exact,
    settings.get('dirichlet_values')
)

# 4. Assemble the global stiffness matrix and body-force vector.
K, F = fem2d.get_global_stiffness_and_force(
    total_nodes, coords, total_element, elements, D_mat, gauss_points,
    gauss_weights, body_force
)

# 5. Add the prescribed boundary traction to F.
F = fem2d.apply_traction_boundary_conditions(
    F, coords, boundary_connect, nodeset_list, settings['traction'],
    settings.get('traction_values'), gauss_points, gauss_weights, exact
)

# 6. Apply displacement constraints and solve KU = F.
K, F = fem2d.apply_essential_boundary_conditions(K, F, fix_nodes)
U = np.linalg.solve(K, F)
print('Problem:', active_problem)
print('Mesh:', settings['mesh_file'])
print('Number of nodes:', total_nodes)
print('Number of Q4 elements:', total_element)
print('Maximum absolute displacement:', np.max(np.abs(U)))

# 7. Save nodal displacements and visualize the result.
np.savetxt(
    output_directory / 'nodal_displacements.txt',
    np.column_stack((coords, U.reshape(-1, 2))),
    header='x  y  u_x  u_y',
    fmt='%.8e',
)
fem2d.plot_mesh(
    coords, elements, settings['title'] + ' Mesh',
    output_directory / 'mesh.png'
)
U_exact = None
exact_stress = None
if has_analytical:
    U_exact = np.concatenate([exact.displacement(x, y) for x, y in coords])
    stress_exact = np.array([exact.stress(x, y) for x, y in coords])
    exact_stress = (stress_exact[:, 0, 0], stress_exact[:, 1, 1], stress_exact[:, 0, 1])

fem2d.plot_displacement_field(
    coords, elements, U, output_directory / 'displacement.png', U_exact
)
sigma11, sigma22, sigma12 = fem2d.compute_stress_at_nodes(
    coords, elements, U, D_mat
)
np.savetxt(
    output_directory / 'nodal_stresses.txt',
    np.column_stack((coords, sigma11, sigma22, sigma12)),
    header='x  y  sigma_11  sigma_22  sigma_12',
    fmt='%.8e',
)
fem2d.plot_stress_contours(
    coords, elements, sigma11, sigma22, sigma12,
    output_directory / 'stress.png', exact_stress
)

# 8. Repeat the same solve for each refinement level when requested.

run_convergence_study = settings.get('convergence', False) and has_analytical
if run_convergence_study:
    fem2d.calculate_convergence(
        settings, input_directory, D_mat, gauss_points, gauss_weights,
        body_force, exact, number_gauss_points, output_directory / 'convergence.png'
    )

print('Results saved to:', output_directory)
