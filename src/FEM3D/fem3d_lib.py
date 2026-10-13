"""Educational H8 3D linear elasticity; SciPy sparse assembly and direct solution.

Voigt order [xx,yy,zz,yz,xz,xy]; strain uses engineering shear.
No large-strain, plasticity, contact, reduced integration, or dynamics.
"""
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
# Outward-oriented face nodes, standard H8 order.
FACES={'XMIN':[0,4,7,3],'XMAX':[1,2,6,5],
       'YMIN':[0,1,5,4],'YMAX':[3,7,6,2],
       'ZMIN':[0,3,2,1],'ZMAX':[4,5,6,7]}
EDGES=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]


def get_elasticity_matrix(E, nu):
    """Build the three-dimensional constitutive matrix."""
    shear_modulus = E / (2.0 * (1.0 + nu))
    lame_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    normal_diagonal = lame_lambda + 2.0 * shear_modulus
    normal_coupling = lame_lambda

    D = np.array([
        [normal_diagonal, normal_coupling, normal_coupling, 0.0, 0.0, 0.0],
        [normal_coupling, normal_diagonal, normal_coupling, 0.0, 0.0, 0.0],
        [normal_coupling, normal_coupling, normal_diagonal, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, shear_modulus, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, shear_modulus, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, shear_modulus],
    ])
    return D


def quadrature(order):
    if not isinstance(order,int) or order<2:raise ValueError('H8 requires integration order >=2')
    return np.polynomial.legendre.leggauss(order)


def shapefunctions_h8(xi, eta, zeta):
    """Return shape functions and derivatives with respect to xi, eta, zeta."""
    N = np.zeros(8)
    dN_dxi = np.zeros(8)
    dN_deta = np.zeros(8)
    dN_dzeta = np.zeros(8)

    for i in range(8):
        xi_i = SIGNS[i, 0]
        eta_i = SIGNS[i, 1]
        zeta_i = SIGNS[i, 2]

        xi_factor = 1.0 + xi_i * xi
        eta_factor = 1.0 + eta_i * eta
        zeta_factor = 1.0 + zeta_i * zeta

        N[i] = xi_factor * eta_factor * zeta_factor / 8.0
        dN_dxi[i] = xi_i * eta_factor * zeta_factor / 8.0
        dN_deta[i] = eta_i * xi_factor * zeta_factor / 8.0
        dN_dzeta[i] = zeta_i * xi_factor * eta_factor / 8.0

    parametric_derivatives = np.array([dN_dxi, dN_deta, dN_dzeta])
    return N, parametric_derivatives


def dNdx(xy, xi, eta, zeta):
    """Calculate physical derivatives and the H8 strain-displacement matrix."""
    N, parametric_derivatives = shapefunctions_h8(xi, eta, zeta)

    # Rows correspond to xi, eta, zeta; columns correspond to x, y, z.
    J = np.zeros((3, 3))
    for i in range(8):
        dN_dxi = parametric_derivatives[0, i]
        dN_deta = parametric_derivatives[1, i]
        dN_dzeta = parametric_derivatives[2, i]
        x_i = xy[i, 0]
        y_i = xy[i, 1]
        z_i = xy[i, 2]

        J[0, 0] += dN_dxi * x_i
        J[0, 1] += dN_dxi * y_i
        J[0, 2] += dN_dxi * z_i
        J[1, 0] += dN_deta * x_i
        J[1, 1] += dN_deta * y_i
        J[1, 2] += dN_deta * z_i
        J[2, 0] += dN_dzeta * x_i
        J[2, 1] += dN_dzeta * y_i
        J[2, 2] += dN_dzeta * z_i

    determinant = np.linalg.det(J)
    if determinant <= 0.0:
        raise ValueError('Nonpositive H8 Jacobian: check element connectivity')

    # Chain rule: J * physical_derivatives = parametric_derivatives.
    physical_derivatives = np.linalg.solve(J, parametric_derivatives)

    B = np.zeros((6, 24))
    for i in range(8):
        dN_dx = physical_derivatives[0, i]
        dN_dy = physical_derivatives[1, i]
        dN_dz = physical_derivatives[2, i]

        ux_column = 3 * i
        uy_column = 3 * i + 1
        uz_column = 3 * i + 2

        # Normal strains: epsilon_xx, epsilon_yy, epsilon_zz.
        B[0, ux_column] = dN_dx
        B[1, uy_column] = dN_dy
        B[2, uz_column] = dN_dz

        # Engineering shear strain gamma_yz = duy/dz + duz/dy.
        B[3, uy_column] = dN_dz
        B[3, uz_column] = dN_dy

        # Engineering shear strain gamma_xz = dux/dz + duz/dx.
        B[4, ux_column] = dN_dz
        B[4, uz_column] = dN_dx

        # Engineering shear strain gamma_xy = dux/dy + duy/dx.
        B[5, ux_column] = dN_dy
        B[5, uy_column] = dN_dx

    return N, B, determinant


def read_lsdyna_k(path):
    """Comma-separated *NODE, *ELEMENT_SOLID (8 nodes), surface *ELEMENT_SHELL.

    Named *SET_NODE_LIST and *SET_SHELL_LIST are merged by the $ comment
    immediately following the keyword. Tetrahedra/degenerate solids unsupported.
    """
    nodes=[];solids=[];faces=[];sets={};section=None;name=None;current=None;entry=None
    for raw in Path(path).read_text().splitlines():
        line=raw.strip()
        if not line:continue
        if line.startswith('*'):
            section=line[1:].upper();name=None;current=None;continue
        if line.startswith('$'):
            if current is None:name=line.lstrip('$#').strip()
            continue
        values=[v.strip() for v in line.split(',') if v.strip()]
        if section=='NODE':nodes.append((int(values[0]),[float(v) for v in values[1:4]]))
        elif section=='ELEMENT_SOLID':
            conn=[int(v) for v in values[2:10]]
            if len(conn)!=8 or len(set(conn))!=8:raise ValueError('Only nondegenerate H8 solids supported')
            solids.append((int(values[0]),conn))
        elif section=='ELEMENT_SHELL':faces.append((int(values[0]),[int(v) for v in values[2:6]]))
        elif section in ('SET_NODE_LIST','SET_SHELL_LIST'):
            if current is None:
                key=name or str(int(values[0]));current=sets.setdefault(key,{'nodes':[],'faces':[]})
                entry='nodes' if section=='SET_NODE_LIST' else 'faces'
            else:current[entry].extend(int(v) for v in values if int(v)!=0)
    node_index={id:i for i,(id,_) in enumerate(nodes)};face_index={id:i for i,(id,_) in enumerate(faces)}
    xy=np.array([v for _,v in nodes],float)
    elements=np.array([[node_index[n] for n in c] for _,c in solids],int).reshape(-1,8)
    surface=np.array([[node_index[n] for n in c] for _,c in faces],int).reshape(-1,4)
    for data in sets.values():
        data['nodes']=np.array([node_index[n] for n in data['nodes']],int)
        data['faces']=np.array([face_index[n] for n in data['faces']],int)
    if not len(elements):raise ValueError('No H8 elements found')
    return xy,elements,surface,sets


def get_element_k_and_f(coord, D_mat, gauss_points, gauss_weights, body_force):
    """Calculate the H8 element stiffness matrix and body-force vector."""
    ke = np.zeros((24, 24))
    fe = np.zeros(24)

    for xi, weight_xi in zip(gauss_points, gauss_weights):
        for eta, weight_eta in zip(gauss_points, gauss_weights):
            for zeta, weight_zeta in zip(gauss_points, gauss_weights):

                # Shape functions, strain-displacement matrix and Jacobian.
                N, B, detJ = dNdx(coord, xi, eta, zeta)

                # Volume integration weight at this Gauss point.
                integration_weight = detJ * weight_xi * weight_eta * weight_zeta

                # Add the Gauss-point contribution to the element stiffness.
                stiffness_contribution = B.T @ D_mat @ B
                ke += stiffness_contribution * integration_weight

                # Map the Gauss point to physical coordinates.
                x_gauss = 0.0
                y_gauss = 0.0
                z_gauss = 0.0
                for i in range(8):
                    x_gauss += N[i] * coord[i, 0]
                    y_gauss += N[i] * coord[i, 1]
                    z_gauss += N[i] * coord[i, 2]

                # Evaluate the three body-force components.
                body_force_vector = body_force(x_gauss, y_gauss, z_gauss)

                # Add the contribution to each node's x, y and z forces.
                for i in range(8):
                    fe[3*i] += N[i] * body_force_vector[0] * integration_weight
                    fe[3*i + 1] += N[i] * body_force_vector[1] * integration_weight
                    fe[3*i + 2] += N[i] * body_force_vector[2] * integration_weight

    return ke, fe


def get_global_stiffness_and_force(xy, elements, D, points, weights, body_force):
    """Assemble H8 contributions in COO format, then convert to CSR."""
    number_nodes = len(xy)
    dofs_per_node = 3
    number_dofs = dofs_per_node * number_nodes
    dofs_per_element = 24

    # Numeric arrays avoid the overhead of millions of Python list objects.
    number_entries = len(elements) * dofs_per_element * dofs_per_element
    row_indices = np.empty(number_entries, dtype=np.int32)
    column_indices = np.empty(number_entries, dtype=np.int32)
    stiffness_values = np.empty(number_entries, dtype=np.float64)
    entry_index = 0
    F = np.zeros(number_dofs)

    for element in elements:
        element_coordinates = xy[element]
        ke, fe = get_element_k_and_f(
            element_coordinates, D, points, weights, body_force
        )

        for i in range(8):
            global_node_i = element[i]
            for component_i in range(dofs_per_node):
                global_row = dofs_per_node * global_node_i + component_i
                local_row = dofs_per_node * i + component_i
                F[global_row] += fe[local_row]

                for j in range(8):
                    global_node_j = element[j]
                    for component_j in range(dofs_per_node):
                        global_column = dofs_per_node * global_node_j + component_j
                        local_column = dofs_per_node * j + component_j

                        row_indices[entry_index] = global_row
                        column_indices[entry_index] = global_column
                        stiffness_values[entry_index] = ke[local_row, local_column]
                        entry_index += 1

    K_coo = coo_matrix(
        (stiffness_values, (row_indices, column_indices)),
        shape=(number_dofs, number_dofs)
    )
    # Shared-node contributions are summed during conversion.
    K = K_coo.tocsr()
    K.eliminate_zeros()
    return K, F


def apply_traction_boundary_conditions(
    F, xy, surface, sets, flags, values, points, weights, exact=None
):
    """Integrate prescribed tractions over outward-oriented Q4 boundary faces."""
    if values is None:
        values = {}

    for name, components in flags.items():
        prescribed = values.get(name, (0.0, 0.0, 0.0))
        face_indices = sets[name]['faces']

        for face_index in face_indices:
            face_nodes = surface[face_index]
            face_coordinates = xy[face_nodes]

            for xi, weight_xi in zip(points, weights):
                for eta, weight_eta in zip(points, weights):
                    N = np.array([
                        (1.0 - xi) * (1.0 - eta) / 4.0,
                        (1.0 + xi) * (1.0 - eta) / 4.0,
                        (1.0 + xi) * (1.0 + eta) / 4.0,
                        (1.0 - xi) * (1.0 + eta) / 4.0,
                    ])
                    dN_dxi = np.array([
                        -(1.0 - eta) / 4.0,
                        (1.0 - eta) / 4.0,
                        (1.0 + eta) / 4.0,
                        -(1.0 + eta) / 4.0,
                    ])
                    dN_deta = np.array([
                        -(1.0 - xi) / 4.0,
                        -(1.0 + xi) / 4.0,
                        (1.0 + xi) / 4.0,
                        (1.0 - xi) / 4.0,
                    ])

                    # Surface tangents, Jacobian and outward normal.
                    tangent_xi = dN_dxi @ face_coordinates
                    tangent_eta = dN_deta @ face_coordinates
                    area_vector = np.cross(tangent_xi, tangent_eta)
                    surface_jacobian = np.linalg.norm(area_vector)
                    normal = area_vector / surface_jacobian
                    integration_weight = surface_jacobian * weight_xi * weight_eta

                    gauss_coordinates = N @ face_coordinates
                    x_gauss = gauss_coordinates[0]
                    y_gauss = gauss_coordinates[1]
                    z_gauss = gauss_coordinates[2]

                    traction = np.zeros(3)
                    if 99 in components:
                        analytical_traction = exact.traction(
                            x_gauss, y_gauss, z_gauss, normal
                        )

                    for component in range(3):
                        flag = components[component]
                        if flag == 1:
                            traction[component] = prescribed[component]
                        elif flag == 99:
                            traction[component] = analytical_traction[component]

                    for i in range(4):
                        global_node = face_nodes[i]
                        for component in range(3):
                            global_dof = 3 * global_node + component
                            nodal_force = N[i] * traction[component]
                            F[global_dof] += nodal_force * integration_weight

    return F


def get_fixnodes_from_sets(sets, flags, xy, exact=None, values=None):
    """Collect prescribed DOF values; use consistent values on shared nodes."""
    if values is None:
        values = {}

    constraints = {}
    for name, components in flags.items():
        prescribed = values.get(name, (0.0, 0.0, 0.0))

        for node in sets[name]['nodes']:
            if 99 in components:
                x_node = xy[node, 0]
                y_node = xy[node, 1]
                z_node = xy[node, 2]
                analytical_displacement = exact.displacement(x_node, y_node, z_node)

            for component in range(3):
                flag = components[component]
                global_dof = 3 * int(node) + component

                if flag == 1:
                    constraints[global_dof] = float(prescribed[component])
                elif flag == 99:
                    constraints[global_dof] = float(analytical_displacement[component])

    return constraints


def solve_constrained(K, F, constraints):
    """Partition the system into free and prescribed displacement DOFs."""
    fixed_dofs = np.array(sorted(constraints), dtype=int)
    all_dofs = np.arange(len(F))
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    U = np.zeros(len(F))
    for dof in fixed_dofs:
        U[dof] = constraints[dof]

    K_free_free = K[free_dofs, :][:, free_dofs]
    K_free_fixed = K[free_dofs, :][:, fixed_dofs]
    prescribed_displacements = U[fixed_dofs]
    effective_force = F[free_dofs] - K_free_fixed @ prescribed_displacements

    if len(free_dofs) > 0:
        U[free_dofs] = spsolve(K_free_free, effective_force)
    reactions = K @ U - F

    return U, reactions, free_dofs


def compute_stress_at_gauss_points(xy, elements, U, D, points):
    """Calculate stress at each H8 Gauss point, in xx,yy,zz,yz,xz,xy order."""
    locations = []
    stresses = []

    for element in elements:
        element_coordinates = xy[element]
        element_displacements = np.zeros(24)
        for i in range(8):
            global_node = element[i]
            for component in range(3):
                local_dof = 3 * i + component
                global_dof = 3 * global_node + component
                element_displacements[local_dof] = U[global_dof]

        for xi in points:
            for eta in points:
                for zeta in points:
                    N, B, detJ = dNdx(element_coordinates, xi, eta, zeta)
                    gauss_coordinates = N @ element_coordinates
                    strain = B @ element_displacements
                    stress = D @ strain
                    locations.append(gauss_coordinates)
                    stresses.append(stress)

    return np.array(locations), np.array(stresses)


def compute_stress_at_nodes(xy, elements, U, D):
    """Average element-corner stresses at each shared node."""
    nodal_stress = np.zeros((len(xy), 6))
    element_count = np.zeros(len(xy))

    for element in elements:
        element_coordinates = xy[element]
        element_displacements = np.zeros(24)
        for i in range(8):
            global_node = element[i]
            for component in range(3):
                local_dof = 3 * i + component
                global_dof = 3 * global_node + component
                element_displacements[local_dof] = U[global_dof]

        for i in range(8):
            xi = SIGNS[i, 0]
            eta = SIGNS[i, 1]
            zeta = SIGNS[i, 2]
            N, B, detJ = dNdx(element_coordinates, xi, eta, zeta)
            strain = B @ element_displacements
            stress = D @ strain
            global_node = element[i]
            nodal_stress[global_node] += stress
            element_count[global_node] += 1

    for node in range(len(xy)):
        nodal_stress[node] /= element_count[node]

    return nodal_stress


def calculate_error_norm(xy, elements, U, D, exact, order=3):
    """Integrate relative displacement L2 and energy errors."""
    points, weights = quadrature(order)
    displacement_error_squared = 0.0
    displacement_exact_squared = 0.0
    energy_error_squared = 0.0
    energy_exact_squared = 0.0
    volume = 0.0

    for element in elements:
        element_coordinates = xy[element]
        element_displacements = np.zeros(24)
        for i in range(8):
            global_node = element[i]
            for component in range(3):
                local_dof = 3 * i + component
                global_dof = 3 * global_node + component
                element_displacements[local_dof] = U[global_dof]
        nodal_displacements = element_displacements.reshape(8, 3)

        for xi, weight_xi in zip(points, weights):
            for eta, weight_eta in zip(points, weights):
                for zeta, weight_zeta in zip(points, weights):
                    N, B, detJ = dNdx(element_coordinates, xi, eta, zeta)
                    integration_weight = detJ * weight_xi * weight_eta * weight_zeta
                    volume += integration_weight

                    gauss_coordinates = N @ element_coordinates
                    x_gauss = gauss_coordinates[0]
                    y_gauss = gauss_coordinates[1]
                    z_gauss = gauss_coordinates[2]

                    exact_displacement = exact.displacement(x_gauss, y_gauss, z_gauss)
                    fem_displacement = N @ nodal_displacements
                    displacement_difference = fem_displacement - exact_displacement

                    exact_tensor = exact.stress(x_gauss, y_gauss, z_gauss)
                    exact_stress = np.array([
                        exact_tensor[0, 0],
                        exact_tensor[1, 1],
                        exact_tensor[2, 2],
                        exact_tensor[1, 2],
                        exact_tensor[0, 2],
                        exact_tensor[0, 1],
                    ])
                    exact_strain = np.linalg.solve(D, exact_stress)
                    fem_strain = B @ element_displacements
                    strain_difference = fem_strain - exact_strain

                    displacement_contribution = displacement_difference @ displacement_difference
                    displacement_reference = exact_displacement @ exact_displacement
                    energy_contribution = strain_difference @ D @ strain_difference
                    energy_reference = exact_strain @ D @ exact_strain

                    displacement_error_squared += displacement_contribution * integration_weight
                    displacement_exact_squared += displacement_reference * integration_weight
                    energy_error_squared += energy_contribution * integration_weight
                    energy_exact_squared += energy_reference * integration_weight

    relative_displacement_error = np.nan
    relative_energy_error = np.nan
    if displacement_exact_squared > 0.0:
        relative_displacement_error = np.sqrt(
            displacement_error_squared / displacement_exact_squared
        )
    if energy_exact_squared > 0.0:
        relative_energy_error = np.sqrt(energy_error_squared / energy_exact_squared)

    relative_errors = np.array([relative_displacement_error, relative_energy_error])
    average_element_volume = volume / len(elements)
    mesh_size = average_element_volume ** (1.0 / 3.0)

    return relative_errors, mesh_size


def solve_problem(settings, input_directory):
    """Repeat the same solution workflow for a selected mesh."""
    mesh_path = Path(input_directory) / settings['mesh_file']
    coords, elements, boundary_faces, sets = read_lsdyna_k(mesh_path)
    E, nu = settings['material_properties']
    D = get_elasticity_matrix(E, nu)

    if 'solution' in settings:
        solution_class = settings['solution']
        parameters = settings['model_characteristics']
        exact = solution_class(*parameters, E, nu)
        body_force = exact.body_force
    else:
        exact = None
        body_force_value = np.array(settings.get('body_force', (0.0, 0.0, 0.0)))

        def body_force(x, y, z):
            return body_force_value

    integration_order = settings.get('integration_order', 2)
    points, weights = quadrature(integration_order)
    K, F = get_global_stiffness_and_force(
        coords, elements, D, points, weights, body_force
    )
    F = apply_traction_boundary_conditions(
        F, coords, boundary_faces, sets,
        settings.get('traction', {}), settings.get('traction_values'),
        points, weights, exact
    )
    constraints = get_fixnodes_from_sets(
        sets, settings['dirichlet'], coords, exact, settings.get('dirichlet_values')
    )
    U, reactions, free_dofs = solve_constrained(K, F, constraints)
    residual_norm = np.linalg.norm(reactions[free_dofs])

    return {
        'coords': coords,
        'elements': elements,
        'U': U,
        'reactions': reactions,
        'F': F,
        'D': D,
        'exact': exact,
        'points': points,
        'residual': float(residual_norm),
        'constraints': constraints,
    }


def write_vtu(path, xy, elements, U, stress):
    """Write readable ASCII VTU: one node or cell per data line."""
    def data_array(name, data, kind='Float64', components=1):
        values = np.asarray(data)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        lines = [
            f'        <DataArray type="{kind}" Name="{name}" '
            f'NumberOfComponents="{components}" format="ascii">'
        ]
        for row in values:
            lines.append('          ' + ' '.join(str(value) for value in row))
        lines.append('        </DataArray>')
        return '\n'.join(lines)

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
        '  <UnstructuredGrid>',
        f'    <Piece NumberOfPoints="{len(xy)}" NumberOfCells="{len(elements)}">',
        '      <Points>',
        data_array('Points', xy, components=3),
        '      </Points>',
        '      <Cells>',
        data_array('connectivity', elements, kind='Int32'),
        data_array('offsets', np.arange(1, len(elements) + 1) * 8, kind='Int32'),
        data_array('types', np.full(len(elements), 12), kind='UInt8'),
        '      </Cells>',
        '      <PointData Vectors="displacement">',
        data_array('displacement', U.reshape(-1, 3), components=3),
        data_array('stress_voigt', stress, components=6),
        '      </PointData>',
        '    </Piece>',
        '  </UnstructuredGrid>',
        '</VTKFile>',
    ]
    Path(path).write_text('\n'.join(lines) + '\n')


def plot_results(path, xy, elements, U, title):
    """Plot original and scaled deformed H8 edges."""
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='3d')
    displacement = U.reshape(-1, 3)
    dimensions = np.ptp(xy, axis=0)
    span = dimensions.max()
    displacement_norms = np.linalg.norm(displacement, axis=1)
    maximum_displacement = displacement_norms.max()
    scale = 1.0
    if maximum_displacement > 0.0:
        scale = 0.08 * span / maximum_displacement

    for element in elements:
        for i, j in EDGES:
            edge_nodes = element[[i, j]]
            original = xy[edge_nodes]
            deformed = original + scale * displacement[edge_nodes]
            ax.plot(original[:, 0], original[:, 1], original[:, 2], color='.7', lw=.3)
            ax.plot(deformed[:, 0], deformed[:, 1], deformed[:, 2], color='tab:blue', lw=.6)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(f'{title}\nDeformation scale = {scale:.3g}')
    ax.set_box_aspect(np.maximum(dimensions, 1e-9))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
