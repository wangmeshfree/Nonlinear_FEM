"""MAE 5036: 3D H8 finite element analysis, step-by-step teaching driver."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import fem3d_lib as fem3d
from fem3d_input_decks import problem_settings


# 1. Select the problem and read its settings.
# active_problem = 'Manufactured_cube'
active_problem = 'Uniaxial_tension'
# active_problem = 'Cantilever_solid'

settings = problem_settings[active_problem]
E, nu = settings['material_properties']
number_gauss_points = settings['integration_order']
has_analytical = 'solution' in settings

if has_analytical:
    exact = settings['solution'](*settings['model_characteristics'], E, nu)
    body_force = exact.body_force
else:
    exact = None
    body_force_value = np.array(settings.get('body_force', (0.0, 0.0, 0.0)))

    def body_force(x, y, z):
        return body_force_value

input_directory = Path(__file__).resolve().parent / 'Inputs'
output_directory = (
    Path(__file__).resolve().parent
    / f"Results_3d_{active_problem}_{Path(settings['mesh_file']).stem}"
)
output_directory.mkdir(exist_ok=True)


# 2. Build the three-dimensional constitutive matrix and Gauss rule.
D_mat = fem3d.get_elasticity_matrix(E, nu)
gauss_points, gauss_weights = fem3d.quadrature(number_gauss_points)


# 3. Read the H8 mesh, boundary faces and named sets.
coords, elements, boundary_faces, nodeset_list = fem3d.read_lsdyna_k(
    input_directory / settings['mesh_file']
)
fix_nodes = fem3d.get_fixnodes_from_sets(
    nodeset_list, settings['dirichlet'], coords, exact,
    settings.get('dirichlet_values')
)


# 4. Assemble the global stiffness matrix and body-force vector.
K, F = fem3d.get_global_stiffness_and_force(
    coords, elements, D_mat, gauss_points, gauss_weights, body_force
)


# 5. Add the prescribed surface traction to F.
F = fem3d.apply_traction_boundary_conditions(
    F, coords, boundary_faces, nodeset_list, settings.get('traction', {}),
    settings.get('traction_values'), gauss_points, gauss_weights, exact
)


# 6. Enforce displacement constraints and solve for the free degrees of freedom.
# Original K and F are retained to calculate reactions.
U, reactions, free_dofs = fem3d.solve_constrained(K, F, fix_nodes)
free_residual = np.linalg.norm(reactions[free_dofs])
print('Problem:', active_problem)
print('Number of nodes:', len(coords))
print('Number of H8 elements:', len(elements))
print('Free-node equilibrium residual:', free_residual)


# 7. Calculate stresses and save numerical results.
nodal_stress = fem3d.compute_stress_at_nodes(coords, elements, U, D_mat)
fem3d.write_vtu(output_directory / 'solution.vtu', coords, elements, U, nodal_stress)


report = {
    'problem': active_problem,
    'nodes': len(coords),
    'elements': len(elements),
    'free_residual': float(free_residual),
    'max_displacement': float(np.max(np.linalg.norm(U.reshape(-1, 3), axis=1))),
    'applied_resultant': F.reshape(-1, 3).sum(axis=0).tolist(),
    'reaction_resultant': reactions.reshape(-1, 3).sum(axis=0).tolist(),
}


# 8. Compare with the analytical solution when available.
if has_analytical:
    errors, h = fem3d.calculate_error_norm(coords, elements, U, D_mat, exact)
    report['relative_displacement_L2'] = float(errors[0])
    report['relative_energy'] = float(errors[1])
    print('Relative displacement L2 error:', errors[0])
    print('Relative energy error:', errors[1])


# 9. Repeat the analysis for the configured refinement levels.
run_convergence_study = settings.get('convergence', False) and has_analytical
if run_convergence_study:
    convergence_rows = []
    for mesh_file in settings['mesh_files']:
        mesh_settings = dict(settings, mesh_file=mesh_file)
        mesh_result = fem3d.solve_problem(mesh_settings, input_directory)
        mesh_errors, h = fem3d.calculate_error_norm(
            mesh_result['coords'], mesh_result['elements'], mesh_result['U'],
            mesh_result['D'], mesh_result['exact']
        )
        convergence_rows.append([
            len(mesh_result['coords']), len(mesh_result['elements']), h,
            mesh_errors[0], mesh_errors[1]
        ])
        print(mesh_file, 'relative L2/energy:', mesh_errors)

    convergence_rows = np.array(convergence_rows)
    np.savetxt(
        output_directory / 'convergence.txt', convergence_rows,
        header='nodes elements h relative_displacement_L2 relative_energy'
    )
    h_values = convergence_rows[:, 2]
    l2_errors = convergence_rows[:, 3]
    energy_errors = convergence_rows[:, 4]
    rates = np.diff(np.log(convergence_rows[:, 3:5]), axis=0) / np.diff(np.log(h_values))[:, None]
    report['convergence_rates_L2_energy'] = rates.tolist()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(h_values, l2_errors, 'o-', label='Relative displacement L2')
    ax.loglog(h_values, energy_errors, 's-', label='Relative energy')
    ax.loglog(h_values, l2_errors[0] * (h_values / h_values[0])**2, '--', label='O(h^2)')
    ax.loglog(h_values, energy_errors[0] * h_values / h_values[0], ':', label='O(h)')
    ax.set(xlabel='h=(volume/elements)^(1/3)', ylabel='Relative error')
    ax.grid()
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_directory / 'convergence.png', dpi=160)
    plt.close(fig)

