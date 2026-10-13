import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def get_elasticity_matrix(E, nu, plane_problem_type):
    """Return the constitutive matrix for plane stress or plane strain."""
    if plane_problem_type == 'plane_stress':
        return E / (1.0 - nu**2) * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])
    if plane_problem_type == 'plane_strain':
        mu = E / (2.0 * (1.0 + nu))
        lambda_ = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        return np.array([
            [lambda_ + 2.0 * mu, lambda_, 0.0],
            [lambda_, lambda_ + 2.0 * mu, 0.0],
            [0.0, 0.0, mu],
        ])
    raise ValueError(f'Unknown 2-D problem type: {plane_problem_type}')

# ===== Extracted from FEM2D notebook cell 10 =====
def read_lsdyna_k(filepath):
    """Parse nodes, shell/beam elements, and *SET_{NODE,BEAM,SHELL}_LIST blocks from an LS-DYNA .k file.

    Returns:
        coords: (nnode, 2) array of x, y node coordinates.
        elements_connect: (nelem, 4) array of 0-based Q4 shell connectivity.
        boundary_connect: (nbeam, 2) array of 0-based beam connectivity (boundary edges).
        nodeset_list: dict physical-group name -> {'ids', 'nodes', 'shells', 'beams'}, where
                      'shells' and 'beams' are 0-based row indices into the connectivity arrays.
    """
    node_ids, node_coords = [], []
    shell_ids, shell_conn = [], []
    beam_ids, beam_conn = [], []
    raw_sets = {}

    section = None
    current_set = None
    current_entry = None
    current_name = None

    with open(filepath, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith('*'):
                section = line[1:].strip()
                current_set = None
                current_name = None
                continue

            if line.startswith('$'):
                # Gmsh writes the physical group name as a comment right after the keyword.
                current_name = line.lstrip('$#').strip() or None
                continue

            values = [v.strip() for v in line.split(',') if v.strip() != '']

            if section == 'NODE':
                node_ids.append(int(values[0]))
                node_coords.append([float(values[1]), float(values[2])])

            elif section == 'ELEMENT_SHELL':
                shell_ids.append(int(values[0]))
                shell_conn.append([int(v) for v in values[2:6]]) # grab the 4 node IDs for the shell element

            elif section == 'ELEMENT_BEAM':
                beam_ids.append(int(values[0]))
                beam_conn.append([int(v) for v in values[2:4]]) # grab the 2 node IDs for the beam element

            elif section in ('SET_NODE_LIST', 'SET_BEAM_LIST', 'SET_SHELL_LIST'):
                if current_set is None:
                    # First data line of a set block is the set id; node and element sets
                    # of the same physical group use different ids, so merge them by name.
                    set_id = int(float(values[0]))
                    key = current_name if current_name else set_id
                    current_set = raw_sets.setdefault(
                        key, {'ids': [], 'nodes': [], 'beams': [], 'shells': []}
                    )
                    current_set['ids'].append(set_id)
                    current_entry = {
                        'SET_NODE_LIST': 'nodes',
                        'SET_BEAM_LIST': 'beams',
                        'SET_SHELL_LIST': 'shells',
                    }[section]
                else:
                    current_set[current_entry].extend(int(v) for v in values)

    node_index = {node_id: i for i, node_id in enumerate(node_ids)}
    shell_index = {elem_id: i for i, elem_id in enumerate(shell_ids)}
    beam_index = {elem_id: i for i, elem_id in enumerate(beam_ids)}

    total_nodes = len(node_ids)
    total_element = len(shell_conn)
    coords = np.array(node_coords)
    elements_connect = np.array([[node_index[n] for n in c] for c in shell_conn], dtype=int)
    boundary_connect = np.array([[node_index[n] for n in c] for c in beam_conn], dtype=int)

    nodeset_list = {
        name: {
            'ids': data['ids'],
            'nodes': np.array([node_index[n] for n in data['nodes']], dtype=int),
            'beams': np.array([beam_index[e] for e in data['beams']], dtype=int),
            'shells': np.array([shell_index[e] for e in data['shells']], dtype=int),
        }
        for name, data in raw_sets.items()
    }

    return total_nodes, coords, total_element, elements_connect, boundary_connect, nodeset_list


def check_set_consistency(nodeset_list, elements, beams):
    """Check that the nodes used by each set's elements reproduce the set node list."""
    report = {}
    for name, data in nodeset_list.items():
        used_nodes = set(beams[data['beams']].ravel().tolist())
        used_nodes |= set(elements[data['shells']].ravel().tolist())
        listed_nodes = set(data['nodes'].tolist())
        report[name] = {
            'consistent': used_nodes == listed_nodes,
            'missing_nodes': np.array(sorted(listed_nodes - used_nodes), dtype=int),
            'extra_nodes': np.array(sorted(used_nodes - listed_nodes), dtype=int),
        }
    return report


def get_fixnodes_from_sets(
    nodeset_list, dirichlet_conditions, coords=None, exact=None, dirichlet_values=None
):
    """Build (node, component, value) constraints.

    Flags are 0 (free), 1 (prescribed value, default 0), or 99 (analytical displacement).
    """
    constraints = []
    for name, flags in dirichlet_conditions.items():
        values = (dirichlet_values or {}).get(name, (0.0, 0.0))
        for node in nodeset_list[name]['nodes']:
            for dof_component, flag in enumerate(flags):
                if flag == 0:
                    continue
                if flag == 1:
                    prescribed_value = values[dof_component]
                elif flag == 99:
                    if coords is None or exact is None:
                        raise ValueError(
                            'Analytical displacement constraints require coords and exact.'
                        )
                    x, y = coords[node]
                    prescribed_value = exact.displacement(x, y)[dof_component]
                else:
                    raise ValueError(f'Unknown Dirichlet flag: {flag}')
                constraints.append([node, dof_component, prescribed_value])

    return np.array(constraints, dtype=float).T




def get_mesh_size(coords, elements):
    """Average element size h = sqrt(domain area / number of elements).

    Used instead of a nominal spacing because the Gmsh meshes are unstructured.
    """
    x = coords[elements][:, :, 0]
    y = coords[elements][:, :, 1]
    areas = 0.5 * np.abs(
        np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
    )
    return np.sqrt(areas.sum() / len(elements))


# ===== Extracted from FEM2D notebook cell 14 =====
def get_gauss_integration_points(number_gauss_points_in_one_direction):
    """Return Gauss-Legendre points and weights on the interval [-1, 1]."""
    if not isinstance(number_gauss_points_in_one_direction, int):
        raise TypeError("The number of Gauss points must be an integer.")
    if number_gauss_points_in_one_direction < 1:
        raise ValueError("The number of Gauss points must be at least 1.")

    return np.polynomial.legendre.leggauss(number_gauss_points_in_one_direction)


def shapefunctions_q4(xi, eta):
    """Return Q4 shape functions and their parametric derivatives."""
    N = np.array([
        0.25 * (1 - xi) * (1 - eta),
        0.25 * (1 + xi) * (1 - eta),
        0.25 * (1 + xi) * (1 + eta),
        0.25 * (1 - xi) * (1 + eta)
    ])
    dNdxi = 0.25 * np.array([
        -(1 - eta),
        (1 - eta),
        (1 + eta),
        -(1 + eta)
    ])
    dNdeta = 0.25 * np.array([
        -(1 - xi),
        -(1 + xi),
        (1 + xi),
        (1 - xi)
    ])
    return N, dNdxi, dNdeta

# ===== Extracted from FEM2D notebook cell 16 =====
def get_element_k_and_f(coord, D_mat, gauss_points, gauss_weights, body_force):
    """Return the Q4 element stiffness matrix and body-force vector."""
    ke = np.zeros((8, 8))
    fe = np.zeros(8)

    for xi, weight_xi in zip(gauss_points, gauss_weights):
        for eta, weight_eta in zip(gauss_points, gauss_weights):

            # Shape functions and their parametric derivatives
            N, dNdxi, dNdeta = shapefunctions_q4(xi, eta)

            # Jacobian matrix 
            # [ \partial x/partial \xi             \partial x/partial \eta ]
            # [ \partial y/partial \xi             \partial y/partial \eta ]
            J = np.zeros((2, 2))
            x_gauss = 0
            y_gauss = 0
            for i in range(4):
                J[0, 0] += dNdxi[i] * coord[i, 0]
                J[0, 1] += dNdeta[i] * coord[i, 0]
                J[1, 0] += dNdxi[i] * coord[i, 1]
                J[1, 1] += dNdeta[i] * coord[i, 1]
                x_gauss += N[i] * coord[i, 0]
                y_gauss += N[i] * coord[i, 1]

            detJ = np.linalg.det(J)
            invJ = np.linalg.inv(J)

            # Compute the derivatives of shape functions with respect to physical coordinates
            dNdxy = invJ.T @ np.array([dNdxi, dNdeta])

            B = np.zeros((3, 8))
            N_matrix = np.zeros((2, 8))
            for i in range(4):
                B[0, 2 * i] = dNdxy[0, i]
                B[1, 2 * i + 1] = dNdxy[1, i]
                B[2, 2 * i] = dNdxy[1, i]
                B[2, 2 * i + 1] = dNdxy[0, i]
                N_matrix[0, 2 * i] = N[i]
                N_matrix[1, 2 * i + 1] = N[i]

            ke += B.T @ D_mat @ B * detJ * weight_xi * weight_eta

            # Physical coordinates of the Gauss point to evaluate the body force.
            fe += N_matrix.T @ body_force(x_gauss, y_gauss) * detJ * weight_xi * weight_eta

    return ke, fe

# ===== Extracted from FEM2D notebook cell 18 =====
def get_global_stiffness_and_force(
    total_nodes, coords, total_element, elements, D_mat, gauss_points, gauss_weights, 
    body_force
):
    """Assemble the global stiffness matrix K and body-force vector F."""
    degrees_of_freedom_per_node = 2
    K = np.zeros((degrees_of_freedom_per_node * total_nodes,
                  degrees_of_freedom_per_node * total_nodes))
    F = np.zeros(degrees_of_freedom_per_node * total_nodes)

    for e in range(total_element):
        # node IDs in the e-th element, like 1 3 5 6 
        element = elements[e]

        # node coordinates for the e-th element
        element_coordinates = coords[element]

        # calculate the elment k and f for the e-th element
        ke, fe = get_element_k_and_f(
            element_coordinates, D_mat, gauss_points, gauss_weights, body_force
        )

        for i in range(4):
            
            for dof_a in range(degrees_of_freedom_per_node):

                # global node id
                row = degrees_of_freedom_per_node * element[i] + dof_a

                # local node id
                local_dof_a = degrees_of_freedom_per_node * i + dof_a

                F[row] += fe[local_dof_a]

                for j in range(4):

                    for dof_b in range(degrees_of_freedom_per_node):

                        # global node id for the b-th local node
                        column = degrees_of_freedom_per_node * element[j] + dof_b

                        # local node id for the b-th local node
                        local_dof_b = degrees_of_freedom_per_node * j + dof_b

                        # local node id for the b-th local node
                        K[row, column] += ke[local_dof_a, local_dof_b]

    return K, F


def apply_traction_boundary_conditions(
    F, coords, boundary_connect, nodeset_list, traction_flags, traction_values,
    gauss_points, gauss_weights, exact
):
    """Add the edge-traction contributions of each loaded boundary set to F.

    Flags per set are (t_x, t_y): 0 = none, 1 = value from traction_values
    (default 0), 99 = analytical traction from exact (needs a solution).
    """
    for set_name, flags in traction_flags.items():
        values = (traction_values or {}).get(set_name, (0.0, 0.0))
        edges = boundary_connect[nodeset_list[set_name]['beams']]

        for edge in edges:
            x1, y1 = coords[edge[0]]
            x2, y2 = coords[edge[1]]
            edge_length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            normal = np.array([y2 - y1, -(x2 - x1)])
            normal = normal / edge_length
            edge_jacobian = edge_length / 2.0

            for gauss_point, gauss_weight in zip(gauss_points, gauss_weights):
                N1 = 0.5 * (1 - gauss_point)
                N2 = 0.5 * (1 + gauss_point)

                # Physical coordinates of the boundary Gauss point.
                x_gauss = N1 * x1 + N2 * x2
                y_gauss = N1 * y1 + N2 * y2

                traction = np.zeros(2)
                for component, flag in enumerate(flags):
                    if flag == 1:
                        traction[component] = values[component]
                    elif flag == 99:
                        if exact is None:
                            raise ValueError('Traction flag 99 (analytical) needs a problem with a solution.')
                        traction[component] = exact.traction(x_gauss, y_gauss, normal)[component]
                    elif flag != 0:
                        raise ValueError(f'Unknown traction flag: {flag}')

                F[2 * edge[0]:2 * edge[0] + 2] += N1 * traction * edge_jacobian * gauss_weight
                F[2 * edge[1]:2 * edge[1] + 2] += N2 * traction * edge_jacobian * gauss_weight

    return F


def apply_essential_boundary_conditions(K, F, fixnodes):
    """Apply prescribed nodal displacements to K and F."""
    degrees_of_freedom_per_node = 2

    for constraint_index in range(fixnodes.shape[1]):
        node = int(fixnodes[0, constraint_index])
        dof_component = int(fixnodes[1, constraint_index])
        prescribed_value = fixnodes[2, constraint_index]
        dof_index = degrees_of_freedom_per_node * node + dof_component

        F -= K[:, dof_index] * prescribed_value
        K[dof_index, :] = 0.0
        K[:, dof_index] = 0.0
        K[dof_index, dof_index] = 1.0
        F[dof_index] = prescribed_value

    return K, F

# ===== Extracted from FEM2D notebook cell 22 =====
def plot_mesh(coords, elements, title="Mesh", output_path=None):
    """Plot mesh"""
    plt.figure(figsize=(8, 8))

    # Plot elements
    for elem in elements:
        x = coords[elem[[0, 1, 2, 3, 0]], 0]
        y = coords[elem[[0, 1, 2, 3, 0]], 1]
        plt.plot(x, y, 'b-', linewidth=0.5)

    # Plot nodes
    plt.plot(coords[:, 0], coords[:, 1], 'ko', markersize=3)

    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()


def compute_stress_at_nodes(coords, elements, U, D_mat):
    
    nnode = len(coords)
    sigma11 = np.zeros(nnode)
    sigma22 = np.zeros(nnode)
    sigma12 = np.zeros(nnode)
    node_count = np.zeros(nnode)
    
    # corners of the reference element
    xi_corners = [-1.0, 1.0, 1.0, -1.0]
    eta_corners = [-1.0, -1.0, 1.0, 1.0]
    
    for elem in elements:
        elem_coords = coords[elem]
        elem_disp = np.zeros(8)
        for i in range(4):
            elem_disp[2*i] = U[2*elem[i]]
            elem_disp[2*i+1] = U[2*elem[i]+1]
        
        # compute stress at each corner
        for corner in range(4):
            xi = xi_corners[corner]
            eta = eta_corners[corner]
            
            N, dNdxi, dNdeta = shapefunctions_q4(xi, eta)
            J = np.zeros((2, 2))
            x_gauss = 0
            y_gauss = 0
            for i in range(4):
                J[0, 0] += dNdxi[i] * elem_coords[i, 0]
                J[0, 1] += dNdeta[i] * elem_coords[i, 0]
                J[1, 0] += dNdxi[i] * elem_coords[i, 1]
                J[1, 1] += dNdeta[i] * elem_coords[i, 1]
                x_gauss += N[i] * elem_coords[i, 0]
                y_gauss += N[i] * elem_coords[i, 1]
                
            detJ = np.linalg.det(J)
            invJ = np.linalg.inv(J)
            dNdxy = invJ.T @ np.array([dNdxi, dNdeta])
            # B matrix
            B = np.zeros((3, 8))
            for i in range(4):
                B[0, 2*i] = dNdxy[0, i]
                B[1, 2*i+1] = dNdxy[1, i]
                B[2, 2*i] = dNdxy[1, i]
                B[2, 2*i+1] = dNdxy[0, i]
            
            # compute stress
            strain = B @ elem_disp
            stress = D_mat @ strain
            
            # accumulate to the node
            node_idx = elem[corner]
            sigma11[node_idx] += stress[0]
            sigma22[node_idx] += stress[1]
            sigma12[node_idx] += stress[2]
            node_count[node_idx] += 1
    
    # average
    sigma11 /= node_count
    sigma22 /= node_count
    sigma12 /= node_count
    
    return sigma11, sigma22, sigma12


def compute_stress_at_gauss_points(coords, elements, U, D_mat, gauss_points):
    """Return physical coordinates and stresses at each element Gauss point."""
    gauss_coordinates = []
    gauss_stresses = []
    nodal_displacements = U.reshape(-1, 2)

    for element in elements:
        element_coords = coords[element]
        element_displacements = nodal_displacements[element].reshape(8)

        for xi in gauss_points:
            for eta in gauss_points:
                shape, dshape_dxi, dshape_deta = shapefunctions_q4(xi, eta)
                jacobian = np.zeros((2, 2))
                for local_node in range(4):
                    jacobian[0, 0] += dshape_dxi[local_node] * element_coords[local_node, 0]
                    jacobian[0, 1] += dshape_deta[local_node] * element_coords[local_node, 0]
                    jacobian[1, 0] += dshape_dxi[local_node] * element_coords[local_node, 1]
                    jacobian[1, 1] += dshape_deta[local_node] * element_coords[local_node, 1]

                dshape_dx = np.linalg.inv(jacobian).T @ np.array(
                    [dshape_dxi, dshape_deta]
                )
                B = np.zeros((3, 8))
                for local_node in range(4):
                    B[0, 2 * local_node] = dshape_dx[0, local_node]
                    B[1, 2 * local_node + 1] = dshape_dx[1, local_node]
                    B[2, 2 * local_node] = dshape_dx[1, local_node]
                    B[2, 2 * local_node + 1] = dshape_dx[0, local_node]

                strain = B @ element_displacements
                gauss_coordinates.append(shape @ element_coords)
                gauss_stresses.append(D_mat @ strain)

    return np.asarray(gauss_coordinates), np.asarray(gauss_stresses)


def get_plot_triangulation(coords, elements):
    """Split each Q4 element into two triangles for contour plotting."""
    triangles = []
    for element in elements:
        triangles.append([element[0], element[1], element[2]])
        triangles.append([element[0], element[2], element[3]])
    return mtri.Triangulation(coords[:, 0], coords[:, 1], triangles)


def plot_displacement_field(coords, elements, U, output_path=None, U_exact=None):
    """Plot u_x and u_y at the nodes; U_exact, if given, is shown in a second row."""
    triangulation = get_plot_triangulation(coords, elements)
    rows = [(U.reshape(-1, 2), 'h')]
    if U_exact is not None:
        rows.append((U_exact.reshape(-1, 2), 'exact'))
    fig, axes = plt.subplots(len(rows), 2, figsize=(10, 5 * len(rows)), squeeze=False)

    levels = [None, None]
    for row_axes, (nodal_values, label) in zip(axes, rows):
        for i, (axis, component, name) in enumerate(zip(row_axes, nodal_values.T, ('x', 'y'))):
            if levels[i] is None:
                levels[i] = np.linspace(component.min(), component.max(), 21)
            contour = axis.tricontourf(triangulation, component, levels=levels[i], cmap='jet', extend='both')
            axis.set_title(rf'$u_{name}^{{{label}}}$' if U_exact is not None else rf'$u_{name}$')
            axis.set_aspect('equal')
            axis.set_xlabel('x')
            axis.set_ylabel('y')
            plt.colorbar(contour, ax=axis)

    plt.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_stress_contours(coords, elements, s11, s22, s12, output_path=None, exact_stress=None):
    """Plot stress components; exact_stress = (s11, s22, s12), if given, is shown in a second row."""
    triangulation = get_plot_triangulation(coords, elements)
    rows = [((s11, s22, s12), 'h')]
    if exact_stress is not None:
        rows.append((exact_stress, 'exact'))
    fig, axes = plt.subplots(len(rows), 3, figsize=(15, 5 * len(rows)), squeeze=False)

    levels = [None, None, None]
    for row_axes, (components, label) in zip(axes, rows):
        for i, (axis, component, name) in enumerate(zip(row_axes, components, ('11', '22', '12'))):
            if levels[i] is None:
                levels[i] = np.linspace(component.min(), component.max(), 21)
            contour = axis.tricontourf(triangulation, component, levels=levels[i], cmap='jet', extend='both')
            axis.set_title(rf'$\sigma_{{{name}}}^{{{label}}}$' if exact_stress is not None else rf'$\sigma_{{{name}}}$')
            axis.set_aspect('equal')
            axis.set_xlabel('x')
            axis.set_ylabel('y')
            plt.colorbar(contour, ax=axis)

    plt.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_gauss_point_stresses(gauss_coordinates, gauss_stresses, output_path=None):
    """Plot stress components as colored samples at the element Gauss points."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    stress_components = [
        (0, r'$\sigma_{11}$'),
        (1, r'$\sigma_{22}$'),
        (2, r'$\sigma_{12}$'),
    ]

    for axis, (component, title) in zip(axes, stress_components):
        scatter = axis.scatter(
            gauss_coordinates[:, 0], gauss_coordinates[:, 1],
            c=gauss_stresses[:, component], s=18, marker='s', cmap='jet',
        )
        axis.set_title(title)
        axis.set_aspect('equal')
        axis.set_xlabel('x')
        axis.set_ylabel('y')
        fig.colorbar(scatter, ax=axis)

    fig.suptitle('Stress at Q4 Gauss Points')
    plt.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_convergence(log_h, log_l2, log_energy, output_path=None):
    """Plot log10(error) against log10(h) so line slopes equal convergence rates."""
    fig, axis = plt.subplots(figsize=(8, 6))
    log_h_ref = np.array([log_h[0], log_h[-1]])
    offset = np.log10(0.9)  # keeps the dashed lines slightly below the data
    axis.plot(log_h, log_l2, 'bo-', label='Relative L2 error')
    axis.plot(
        log_h_ref, log_l2[0] + offset + 2.0 * (log_h_ref - log_h[0]),
        'b--', label='Slope 2',
    )
    axis.plot(log_h, log_energy, 'rs-', label='Relative Energy error')
    axis.plot(
        log_h_ref, log_energy[0] + offset + 1.0 * (log_h_ref - log_h[0]),
        'r--', label='Slope 1',
    )
    axis.set_xlabel(r'$\log_{10}(h)$,  $h = \sqrt{A/n_{elem}}$')
    axis.set_ylabel(r'$\log_{10}$(relative error)')
    axis.grid(True, alpha=0.3)
    axis.legend()
    axis.set_title('Q4 Displacement Convergence')
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()



# ===== Extracted from FEM2D notebook cell 26 =====
def calculate_error_norm(coords, elements, U, D_mat, exact, gauss_points):
    
    print("\n===== Calculate error norm =====")
    
    # Four-point quadrature also resolves displacement norms accurately.
    gauss_pts, weights = get_gauss_integration_points(gauss_points)
    
    energy_error = 0.0  # integrated squared energy-norm error
    energy_exact = 0.0
    disp_error = 0.0
    disp_exact = 0.0

    S = np.linalg.inv(D_mat)
    
    # Loop over each element
    for elem in elements:
        elem_coords = coords[elem]
        elem_disp = np.zeros(8)
        for i in range(4):
            elem_disp[2*i] = U[2*elem[i]]
            elem_disp[2*i+1] = U[2*elem[i]+1]
        
        # Integrate over each Gauss point within the element
        for xi in gauss_pts:
            for eta in gauss_pts:
                # Shape functions and their derivatives
                N, dNdxi, dNdeta = shapefunctions_q4(xi, eta)
                
                # Jacobian matrix
                J = np.array([dNdxi, dNdeta]) @ elem_coords
                detJ = np.linalg.det(J)
                invJ = np.linalg.inv(J)
                dNdxy = invJ @ np.array([dNdxi, dNdeta])
                
                
                B = np.zeros((3, 8))
                for i in range(4):
                    B[0, 2*i] = dNdxy[0, i]
                    B[1, 2*i+1] = dNdxy[1, i]
                    B[2, 2*i] = dNdxy[1, i]
                    B[2, 2*i+1] = dNdxy[0, i]
                
                
                strain_fem = B @ elem_disp
                
                
                x_gp = N @ elem_coords[:, 0]
                y_gp = N @ elem_coords[:, 1]
                
                
                stress_cart = exact.stress(x_gp, y_gp)
                stress_exact = np.array([stress_cart[0,0], stress_cart[1,1], stress_cart[0,1]])
                
                
                strain_exact = S @ stress_exact
                
                strain_error = strain_fem - strain_exact
                

                energy_error += strain_error @ D_mat @ strain_error * detJ * weights[0] * weights[1]
                energy_exact += (
                    strain_exact @ D_mat @ strain_exact
                    * detJ * weights[0] * weights[1]
                )



                exact_displacement = exact.displacement(x_gp, y_gp)
                ux_gp = N @ elem_disp.reshape(4, 2)[:, 0]
                uy_gp = N @ elem_disp.reshape(4, 2)[:, 1]

                displacement_error = np.array([ux_gp, uy_gp]) - exact_displacement
                disp_error += displacement_error @ displacement_error * detJ * weights[0] * weights[1]
                disp_exact += exact_displacement @ exact_displacement * detJ * weights[0] * weights[1]

    
    # Calculate the energy norm error (square root)
    energy = np.sqrt(energy_error) / np.sqrt(energy_exact)
    L2 = np.sqrt(disp_error) / np.sqrt(disp_exact) 


    return L2, energy


def calculate_convergence(
    settings, input_directory, D_mat, gauss_points, gauss_weights, body_force,
    exact, number_gauss_points, output_path=None
):
    """Solve each mesh in settings['mesh_files'], print the errors and rates, and plot them."""
    h_vals = []
    l2_error_vals = []
    energy_error_vals = []
    for mesh_file in settings['mesh_files']:
        nnode, mesh_coords, nelem, mesh_elements, mesh_boundary, mesh_sets = (
            read_lsdyna_k(input_directory / mesh_file)
        )
        mesh_fix = get_fixnodes_from_sets(
            mesh_sets, settings['dirichlet'], mesh_coords, exact,
            settings.get('dirichlet_values')
        )
        mesh_K, mesh_F = get_global_stiffness_and_force(
            nnode, mesh_coords, nelem, mesh_elements, D_mat, gauss_points,
            gauss_weights, body_force
        )
        mesh_F = apply_traction_boundary_conditions(
            mesh_F, mesh_coords, mesh_boundary, mesh_sets, settings['traction'],
            settings.get('traction_values'), gauss_points, gauss_weights, exact
        )
        mesh_K, mesh_F = apply_essential_boundary_conditions(
            mesh_K, mesh_F, mesh_fix
        )
        mesh_U = np.linalg.solve(mesh_K, mesh_F)

        h = get_mesh_size(mesh_coords, mesh_elements)

        relative_l2_error, relative_energy_error = calculate_error_norm(
            mesh_coords, mesh_elements, mesh_U, D_mat, exact, number_gauss_points
        )

        h_vals.append(h)
        l2_error_vals.append(relative_l2_error)
        energy_error_vals.append(relative_energy_error)

        print(f'{mesh_file}:')
        print(f'  Nodes: {nnode}; elements: {nelem}')
        print(f'  h = {h:.4f}')
        print(f'  Relative L2 error = {relative_l2_error:.2%}')
        print(f'  Relative Energy error = {relative_energy_error:.2%}')
        print()

    log_h_vals = np.log10(h_vals)
    log_l2_errors = np.log10(l2_error_vals)
    log_energy_errors = np.log10(energy_error_vals)
    l2_rates = np.diff(log_l2_errors) / np.diff(log_h_vals)
    energy_rates = np.diff(log_energy_errors) / np.diff(log_h_vals)
    print(f'Average L2 convergence rate: {np.mean(l2_rates):.4f}')
    print(f'Average Energy convergence rate: {np.mean(energy_rates):.4f}')
    plot_convergence(log_h_vals, log_l2_errors, log_energy_errors, output_path) 
