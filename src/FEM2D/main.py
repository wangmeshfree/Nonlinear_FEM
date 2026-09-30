"""2D Q4 FEM driver: select a problem, solve, plot, and save results."""

from pathlib import Path

import numpy as np

import fem2d_lib as fem2d


# 1. Choose the problem and material.
active_problem = 'Cantilever_beam'
run_convergence_study = True

problem_settings = {
    'Plate_with_hole': {
        'model_characteristics': [1.0, 1.0],  # [hole radius R, tensile load Tx]
        'mesh_file': '2d_plate_with_hole/model_plate_with_hole_level_2.k',
        'mesh_files': [
            '2d_plate_with_hole/model_plate_with_hole_level_1.k',
            '2d_plate_with_hole/model_plate_with_hole_level_2.k',
            '2d_plate_with_hole/model_plate_with_hole_level_3.k',
        ],
        'material_properties': [1.0e3, 0.30],
        'integration_order': 2,
        'dirichlet': {'FIXED_LEFT': (1, 0), 'FIXED_BOTTOM': (0, 1)},
        'traction': ['TRACTION_RIGHT', 'TRACTION_TOP'],
        'title': 'Plate with Hole',
        '2d_problem_type': 'plane_strain',
    },
    'Cantilever_beam': {
        'model_characteristics': [10.0, 2.0, 5.0],  # [length L, height H, tip load P]
        'mesh_file': '2d_cantilever_beam/model_cantilever_level_2.k',
        'mesh_files': [
            '2d_cantilever_beam/model_cantilever_level_1.k',
            '2d_cantilever_beam/model_cantilever_level_2.k',
            '2d_cantilever_beam/model_cantilever_level_3.k',
        ],
        'material_properties': [1.0e3, 0.30],
        'integration_order': 2,
        'dirichlet': {'FIXED_LEFT': (99, 99)},
        'traction': ['TRACTION_RIGHT'],
        'title': 'Cantilever Beam',
        '2d_problem_type': 'plane_stress',
    },
}

settings = problem_settings[active_problem]
model_characteristics = settings['model_characteristics']
E, nu = settings['material_properties']
number_gauss_points_in_one_direction = settings['integration_order']
plane_problem_type = settings['2d_problem_type']
input_directory = Path(__file__).resolve().parent / 'Inputs'
output_directory = (
    Path(__file__).resolve().parent
    / f"Results_2d_{active_problem}_{Path(settings['mesh_file']).stem}"
)
output_directory.mkdir(exist_ok=True)


# 2. Build the constitutive matrix and Gauss rule.
D_mat = fem2d.get_elasticity_matrix(E, nu, plane_problem_type)
gauss_points, gauss_weights = fem2d.get_gauss_integration_points(
    number_gauss_points_in_one_direction
)


# 3. Read the selected LS-DYNA keyword mesh and its named boundary sets.
total_nodes, coords, total_element, elements, boundary_connect, nodeset_list = (
    fem2d.read_lsdyna_k(input_directory / settings['mesh_file'])
)
fix_nodes = fem2d.get_fixnodes_from_sets(
    nodeset_list, settings['dirichlet'], coords, model_characteristics,
    [E, nu], plane_problem_type, active_problem
)
fix_nodes = fem2d.add_cantilever_center_constraint(
    fix_nodes, coords, active_problem, settings['dirichlet']
)
traction_elements = fem2d.get_traction_element_indices(
    nodeset_list, settings['traction']
)


# 4. Assemble the global stiffness matrix and body-force vector.
K, F = fem2d.get_global_stiffness_and_force(
    total_nodes, coords, total_element, elements, D_mat, gauss_points,
    gauss_weights, model_characteristics, active_problem
)


# 5. Add the prescribed boundary traction to F.
F = fem2d.apply_traction_boundary_conditions(
    F, coords, boundary_connect, traction_elements, gauss_points,
    gauss_weights, model_characteristics, active_problem
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
fem2d.plot_displacement_field(
    coords, elements, U, output_directory / 'displacement.png'
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
    output_directory / 'stress.png'
)


# 8. Compare with the analytical reference stress.
energy_error = fem2d.calculate_energy_norm_error(
    coords, elements, U, D_mat, model_characteristics, active_problem
)
print('Main mesh energy norm error:', energy_error)


# 9. Repeat the same solve for each refinement level when requested.
if run_convergence_study:
    h_vals = []
    err_vals = []
    convergence_rows = []
    for mesh_file in settings['mesh_files']:
        nnode, mesh_coords, nelem, mesh_elements, mesh_boundary, mesh_sets = (
            fem2d.read_lsdyna_k(input_directory / mesh_file)
        )
        mesh_fix = fem2d.get_fixnodes_from_sets(
            mesh_sets, settings['dirichlet'], mesh_coords, model_characteristics,
            [E, nu], plane_problem_type, active_problem
        )
        mesh_fix = fem2d.add_cantilever_center_constraint(
            mesh_fix, mesh_coords, active_problem, settings['dirichlet']
        )
        mesh_traction = fem2d.get_traction_element_indices(
            mesh_sets, settings['traction']
        )
        mesh_K, mesh_F = fem2d.get_global_stiffness_and_force(
            nnode, mesh_coords, nelem, mesh_elements, D_mat, gauss_points,
            gauss_weights, model_characteristics, active_problem
        )
        mesh_F = fem2d.apply_traction_boundary_conditions(
            mesh_F, mesh_coords, mesh_boundary, mesh_traction, gauss_points,
            gauss_weights, model_characteristics, active_problem
        )
        mesh_K, mesh_F = fem2d.apply_essential_boundary_conditions(
            mesh_K, mesh_F, mesh_fix
        )
        mesh_U = np.linalg.solve(mesh_K, mesh_F)
        h = fem2d.get_mesh_size(mesh_coords, mesh_elements)
        error = fem2d.calculate_energy_norm_error(
            mesh_coords, mesh_elements, mesh_U, D_mat,
            model_characteristics, active_problem
        )
        h_vals.append(h)
        err_vals.append(error)
        convergence_rows.append((nnode, nelem, h, error))
        print(f'{mesh_file}: {nnode} nodes, {nelem} elements, '
              f'h = {h:.4f}, energy error = {error:.6e}')

    np.savetxt(
        output_directory / 'convergence.txt',
        convergence_rows,
        header='nodes  Q4_elements  h  energy_norm_error',
        fmt=['%d', '%d', '%.8e', '%.8e'],
    )
    for i in range(len(h_vals) - 1):
        rate = np.log(err_vals[i + 1] / err_vals[i]) / np.log(
            h_vals[i + 1] / h_vals[i]
        )
        print(f'Convergence rate {i + 1}: {rate:.4f}')
    fem2d.plot_convergence(
        h_vals, err_vals, output_directory / 'convergence.png'
    )

print('Results saved to:', output_directory)
