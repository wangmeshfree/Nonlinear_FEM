import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# Analytical fields and loads
def get_stress_at_point(x, y, model_characteristics, problem_name):
    """Return the analytical reference stress for the selected problem."""
    if problem_name == 'Plate_with_hole':
        return stress_plate_with_hole(x, y, model_characteristics)
    elif problem_name == 'Cantilever_beam':
        return stress_cantilever_beam(x, y, model_characteristics)
    else:
        raise ValueError(f"Please tell me the analytical stress field for problem: {problem_name}")


def get_traction_at_integration_point(x, y, model_characteristics, normal, problem_name):
    """Return the traction vector at a boundary integration point."""
    
    analytical_problems = ('Plate_with_hole', 'Cantilever_beam')

    if problem_name in analytical_problems:
        stress = get_stress_at_point(x, y, model_characteristics, problem_name)
        return stress @ normal
    else:
        raise ValueError(f"Please specify how the traction is prescribed for problem: {problem_name}")


def get_body_force(x, y, model_characteristics, problem_name):
    if problem_name == 'Plate_with_hole':
        return body_force_plate_with_hole(x, y, model_characteristics)
    elif problem_name == 'Cantilever_beam':
        return body_force_cantilever_beam(x, y, model_characteristics)
    else:
        return np.array([0.0, 0.0])

# ===== Extracted from FEM2D notebook cell 8 =====
#########################################################
# Analytical stress and body-force fields used by the FEM driver.
#
# The cantilever convention is: x in [0, L_b], y in [-H/2, H/2],
# unit out-of-plane thickness, and a positive tip load P in +y.
# The exact field below is the classical Euler--Bernoulli beam field.
#########################################################
def stress_plate_with_hole(x, y, model_characteristics):
    """Return the analytical (Kirsch) stress tensor [[sxx, sxy], [sxy, syy]] at (x, y).

    Points inside the hole are clamped to r = R so the polar formulas stay valid.
    """
    R = model_characteristics[0]
    Tx = model_characteristics[1]
    r = max(np.sqrt(x**2 + y**2), R)
    th = np.arctan2(y, x)
    rr = R / r

    # Polar components of the plate-with-hole solution.
    srr = 0.5*Tx*(1.0 - rr**2.0) + 0.5*Tx*(1.0 - 4.0*rr**2.0 + 3.0*rr**4.0)*np.cos(2.0*th)
    stt = 0.5*Tx*(1.0 + rr**2.0) - 0.5*Tx*(1.0 + 3.0*rr**4.0)*np.cos(2.0*th)
    srt = -0.5*Tx*(1.0 + 2.0*rr**2.0 - 3.0*rr**4.0)*np.sin(2.0*th)

    # Rotate the polar stress tensor into Cartesian components.
    RM = np.array([[np.cos(th), np.sin(th)],
                   [-np.sin(th), np.cos(th)]])
    return RM.T @ np.array([[srr, srt], [srt, stt]]) @ RM
def body_force_plate_with_hole(x, y, model_characteristics):
    """Return the body-force components at the physical point (x, y) for a plate with a hole."""
    bx = 0.0
    by = 0.0
    return np.array([bx, by])
def stress_cantilever_beam(x, y, model_characteristics):
    """Return the classical cantilever stress tensor at (x, y).

    ``model_characteristics = [L_b, H, P]`` with unit thickness.
    The beam occupies 0 <= x <= L_b and -H/2 <= y <= H/2.
    A positive P acts in the +y direction at x = L_b.

    The field is
        sigma_xx = -M(x) * y / I,
        sigma_xy = 3*P/(2*H) * (1 - 4*y^2/H^2),
        sigma_yy = 0,
    where M(x) = P*(L_b-x) and I = H^3/12.
    It satisfies equilibrium and gives zero traction on y = +/- H/2.
    At x = L_b, its shear traction integrates to the tip load P.
    """
    beam_length, beam_height, tip_load = model_characteristics[:3]
    y_centred = y
    area = beam_height          # unit thickness
    second_moment = beam_height**3 / 12.0
    bending_moment = tip_load * (beam_length - x)

    sxx = -bending_moment * y_centred / second_moment
    syy = 0.0
    sxy = (3.0 * tip_load / (2.0 * area)) * (
        1.0 - 4.0 * y_centred**2 / beam_height**2
    )
    return np.array([[sxx, sxy], [sxy, syy]])

def body_force_cantilever_beam(x, y, model_characteristics):
    """Return the body-force components at the physical point (x, y) for a cantilever beam."""
    bx = 0.0
    by = 0.0
    return np.array([bx, by])


def get_displacement_at_point(
    x, y, model_characteristics, model_materials, plane_problem_type, problem_name
):
    """Return the analytical displacement for the selected benchmark."""
    if problem_name == 'Cantilever_beam':
        return displacement_cantilever_beam(
            x, y, model_characteristics, model_materials, plane_problem_type
        )
    raise ValueError(
        f'No analytical displacement field is defined for {problem_name}'
    )


def displacement_cantilever_beam(
    x, y, model_characteristics, model_materials, plane_problem_type
):
    """2-D cantilever displacement solution used for analytical BCs."""
    beam_length, beam_height, tip_load = model_characteristics[:3]
    E, nu = model_materials
    if plane_problem_type == 'plane_stress':
        E_eff, nu_eff = E, nu
    elif plane_problem_type == 'plane_strain':
        E_eff = E / (1.0 - nu**2)
        nu_eff = nu / (1.0 - nu)
    else:
        raise ValueError(f'Unknown plane problem type: {plane_problem_type}')

    second_moment = beam_height**3 / 12.0
    ux = -tip_load * y / (6.0 * E_eff * second_moment) * (
        (6.0 * beam_length - 3.0 * x) * x
        + (2.0 + nu_eff) * (y**2 - beam_height**2 / 4.0)
    )
    uy = tip_load / (6.0 * E_eff * second_moment) * (
        3.0 * nu_eff * y**2 * (beam_length - x)
        + (4.0 + 5.0 * nu_eff) * beam_height**2 * x / 4.0
        + (3.0 * beam_length - x) * x**2
    )
    return np.array([ux, uy])


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
    nodeset_list, dirichlet_conditions, coords=None, model_characteristics=None,
    model_materials=None, plane_problem_type=None, problem_name=None
):
    """Build (node, component, value) constraints.

    Flags are 0 (free), 1 (zero displacement), or 99 (analytical displacement).
    """
    constraints = []
    for name, flags in dirichlet_conditions.items():
        for node in nodeset_list[name]['nodes']:
            for dof_component, flag in enumerate(flags):
                if flag == 0:
                    continue
                if flag == 1:
                    prescribed_value = 0.0
                elif flag == 99:
                    if any(value is None for value in (
                        coords, model_characteristics, model_materials,
                        plane_problem_type, problem_name
                    )):
                        raise ValueError(
                            'Analytical displacement constraints require model data.'
                        )
                    x, y = coords[node]
                    prescribed_value = get_displacement_at_point(
                        x, y, model_characteristics, model_materials,
                        plane_problem_type, problem_name
                    )[dof_component]
                else:
                    raise ValueError(f'Unknown Dirichlet flag: {flag}')
                constraints.append([node, dof_component, prescribed_value])

    return np.array(constraints, dtype=float).T


def add_cantilever_center_constraint(
    fix_nodes, coords, problem_name, dirichlet_conditions
):
    """Add the old centre constraint only for zero-value beam BCs."""
    if problem_name != 'Cantilever_beam':
        return fix_nodes
    if any(flag == 99 for flags in dirichlet_conditions.values() for flag in flags):
        return fix_nodes
    left_nodes = np.where(np.isclose(coords[:, 0], np.min(coords[:, 0])))[0]
    center_node = left_nodes[np.argmin(np.abs(coords[left_nodes, 1]))]
    extra = np.array([[center_node], [1], [0.0]], dtype=float)
    return np.column_stack((fix_nodes, extra))

def get_traction_element_indices(nodeset_list, traction_sets):
    """Return the boundary element row indices belonging to the named traction sets."""
    return np.concatenate([nodeset_list[name]['beams'] for name in traction_sets])

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

# ===== Extracted from FEM2D notebook cell 12 =====
# Quarter plate with a circular hole: x >= 0, y >= 0
# The outer edges receive the analytical traction from the plate-with-hole solution.
def get_mesh_platehole(nelem_1D, R, L):
    nnode_1D = nelem_1D + 1
    nnode = nnode_1D**2 * 2 - nnode_1D
    nelem = nelem_1D**2 * 2
    tolerance = 1.0e-10

    th = np.linspace(0.0, np.pi / 4.0, nnode_1D)
    x_1, y_1, x_2, y_2 = [], [], [], []

    for irow in range(nnode_1D):
        x_row = np.linspace(np.cos(th[irow]) * R, L, nnode_1D)
        y_row = np.linspace(np.sin(th[irow]) * R, L / nelem_1D * irow, nnode_1D)
        x_1.extend(x_row)
        y_1.extend(y_row)

        if irow != nnode_1D - 1:
            x_2 = list(y_row) + x_2
            y_2 = list(x_row) + y_2

    coords = np.array([x_1 + x_2, y_1 + y_2])

    connect = np.zeros((4, nelem), dtype=int)
    rowcount = 0
    for elementcount in range(nelem):
        connect[0, elementcount] = elementcount + rowcount
        connect[1, elementcount] = elementcount + rowcount + 1
        connect[2, elementcount] = elementcount + rowcount + nnode_1D + 1
        connect[3, elementcount] = elementcount + rowcount + nnode_1D

        if (elementcount + 1) % nelem_1D == 0:
            rowcount += 1

    # Symmetry boundary conditions: u_y(x, 0) = 0 and u_x(0, y) = 0.
    nodes_on_x_axis = np.where(np.isclose(coords[1, :], 0.0, atol=tolerance))[0]
    nodes_on_y_axis = np.where(np.isclose(coords[0, :], 0.0, atol=tolerance))[0]
    nfix = len(nodes_on_x_axis) + len(nodes_on_y_axis)
    fixnodes = np.zeros((3, nfix), dtype=int)
    fixnodes[0, :len(nodes_on_x_axis)] = nodes_on_x_axis
    fixnodes[1, :len(nodes_on_x_axis)] = 1
    fixnodes[0, len(nodes_on_x_axis):] = nodes_on_y_axis
    fixnodes[1, len(nodes_on_x_axis):] = 0

    # Q4 local faces, listed counter-clockwise around the reference element.
    q4_faces = {
        0: np.array([0, 1]),
        1: np.array([1, 2]),
        2: np.array([2, 3]),
        3: np.array([3, 0]),
    }

    # Each entry defines a physical traction boundary and its outward normal.
    traction_boundaries = [
        {"coordinate": 0, "value": L, "normal": np.array([1.0, 0.0])},
        {"coordinate": 1, "value": L, "normal": np.array([0.0, 1.0])},
    ]

    loaded_faces = []
    for element_index, element in enumerate(connect.T):
        for face, local_nodes in q4_faces.items():
            face_nodes = element[local_nodes]

            for boundary in traction_boundaries:
                face_coordinates = coords[boundary["coordinate"], face_nodes]
                if np.all(np.isclose(face_coordinates, boundary["value"], atol=tolerance)):
                    loaded_faces.append((element_index, face, boundary["normal"]))

    dloads = np.zeros((4, len(loaded_faces)))
    for load_index, (element_index, face, normal) in enumerate(loaded_faces):
        dloads[0, load_index] = element_index
        dloads[1, load_index] = face
        dloads[2:4, load_index] = normal

    return nnode, coords.T, nelem, connect.T, fixnodes, dloads

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
def get_element_k_and_f(coord, D_mat, gauss_points, gauss_weights, model_characteristics, problem_name):
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
            body_force = get_body_force(x_gauss, y_gauss, model_characteristics, problem_name)
            fe += N_matrix.T @ body_force * detJ * weight_xi * weight_eta

    return ke, fe

# ===== Extracted from FEM2D notebook cell 18 =====
def get_global_stiffness_and_force(
    total_nodes, coords, total_element, elements, D_mat, gauss_points, gauss_weights, 
    model_characteristics, problem_name
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
            element_coordinates, D_mat, gauss_points, gauss_weights, model_characteristics, problem_name
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
    F, coords, boundary_connect, traction_elements, gauss_points, gauss_weights,
    model_characteristics, problem_name
):
    """Add the edge-traction contributions of the loaded boundary elements to F.

    The traction is evaluated at each boundary Gauss point from the outward normal
    of that edge.
    """
    edges = boundary_connect[traction_elements]

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

            traction = get_traction_at_integration_point(
                x_gauss, y_gauss, model_characteristics, normal, problem_name
            )

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
    
    # 参考单元的角点
    xi_corners = [-1.0, 1.0, 1.0, -1.0]
    eta_corners = [-1.0, -1.0, 1.0, 1.0]
    
    for elem in elements:
        elem_coords = coords[elem]
        elem_disp = np.zeros(8)
        for i in range(4):
            elem_disp[2*i] = U[2*elem[i]]
            elem_disp[2*i+1] = U[2*elem[i]+1]
        
        # 计算每个角点的应力
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
            # B矩阵
            B = np.zeros((3, 8))
            for i in range(4):
                B[0, 2*i] = dNdxy[0, i]
                B[1, 2*i+1] = dNdxy[1, i]
                B[2, 2*i] = dNdxy[1, i]
                B[2, 2*i+1] = dNdxy[0, i]
            
            # 计算应力
            strain = B @ elem_disp
            stress = D_mat @ strain
            
            # 累加到节点
            node_idx = elem[corner]
            sigma11[node_idx] += stress[0]
            sigma22[node_idx] += stress[1]
            sigma12[node_idx] += stress[2]
            node_count[node_idx] += 1
    
    # 平均
    sigma11 /= node_count
    sigma22 /= node_count
    sigma12 /= node_count
    
    return sigma11, sigma22, sigma12

def get_plot_triangulation(coords, elements):
    """Split each Q4 element into two triangles for contour plotting."""
    triangles = []
    for element in elements:
        triangles.append([element[0], element[1], element[2]])
        triangles.append([element[0], element[2], element[3]])
    return mtri.Triangulation(coords[:, 0], coords[:, 1], triangles)


def plot_displacement_field(coords, elements, U, output_path=None):
    """Plot the x- and y-displacement components at the mesh nodes."""
    nodal_displacements = U.reshape(-1, 2)
    triangulation = get_plot_triangulation(coords, elements)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    displacement_plots = [
        (nodal_displacements[:, 0], r'$u_x$'),
        (nodal_displacements[:, 1], r'$u_y$')
    ]

    for axis, (displacement_component, title) in zip(axes, displacement_plots):
        contour = axis.tricontourf(
            triangulation, displacement_component, levels=20, cmap='jet'
        )
        axis.set_title(title)
        axis.set_aspect('equal')
        axis.set_xlabel('x')
        axis.set_ylabel('y')
        plt.colorbar(contour, ax=axis)

    plt.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

def plot_stress_contours(coords, elements, s11, s22, s12, output_path=None):
    """Plot stress contours"""
    triangulation = get_plot_triangulation(coords, elements)

    # Plot each stress component
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Sigma_11
    tcf1 = axes[0].tricontourf(triangulation, s11, levels=20, cmap='jet')
    axes[0].set_title(r'$\sigma_{11}$')
    axes[0].set_aspect('equal')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(tcf1, ax=axes[0])

    # Sigma_22
    tcf2 = axes[1].tricontourf(triangulation, s22, levels=20, cmap='jet')
    axes[1].set_title(r'$\sigma_{22}$')
    axes[1].set_aspect('equal')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    plt.colorbar(tcf2, ax=axes[1])

    # Sigma_12
    tcf3 = axes[2].tricontourf(triangulation, s12, levels=20, cmap='jet')
    axes[2].set_title(r'$\sigma_{12}$')
    axes[2].set_aspect('equal')
    axes[2].set_xlabel('x')
    axes[2].set_ylabel('y')
    plt.colorbar(tcf3, ax=axes[2])

    plt.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_convergence(h_vals, err_vals, output_path=None):
    """Plot energy error against mesh size with an O(h) reference line."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(h_vals, err_vals, 'bo-', label='FEM Error')
    reference = err_vals[0] * np.array([h_vals[0], h_vals[-1]]) / h_vals[0]
    ax.loglog([h_vals[0], h_vals[-1]], reference, 'r--', label=r'$O(h)$ Reference')
    ax.set_xlabel(r'Mesh size $h = \sqrt{A/n_{elem}}$')
    ax.set_ylabel('Energy norm error')
    ax.set_title('Q4 Element Convergence')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

# ===== Extracted from FEM2D notebook cell 26 =====
def calculate_energy_norm_error(coords, elements, U, D_mat, model_characteristics, problem_name):
    
    print("\n===== Calculate error norm =====")
    
    # Gauss points and weights
    gauss_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]
    weights = [1.0, 1.0]
    
    enorm2 = 0.0  # Square of the energy norm error
    
    
    S = np.linalg.inv(D_mat)
    # print(f" S = \n{S}\n")
    
    # Loop over each element
    for elem_idx, elem in enumerate(elements):
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
                
                
                stress_cart = get_stress_at_point(x_gp, y_gp, model_characteristics, problem_name)
                stress_exact = np.array([stress_cart[0,0], stress_cart[1,1], stress_cart[0,1]])
                
                
                strain_exact = S @ stress_exact
                
                
                strain_error = strain_fem - strain_exact
                
                
                contribution = strain_error @ D_mat @ strain_error * detJ * weights[0] * weights[1]
                enorm2 += contribution
                
                
                # if elem_idx % 100 == 0 and xi == gauss_pts[0] and eta == gauss_pts[0]:
                #     print(f"Element {elem_idx}, Gauss point ({xi:.2f},{eta:.2f}):")
                #     print(f"  Strain error: {strain_error}")
                #     print(f"  Contribution: {contribution:.6f}, Cumulative energy error squared: {enorm2:.6f}\n")
    
    # Calculate the energy norm error (square root)
    energy_norm_error = np.sqrt(enorm2)
    print(f"===== Error norm calculation completed =====")
    print(f"Total energy norm error = {energy_norm_error:.6f}\n")
    
    return energy_norm_error
