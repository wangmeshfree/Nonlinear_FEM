# ******************************************************************************
# *                                                                            *
# *        MAE 5036 - Advanced Computational Solid Mechanics II                *
# *                                                                            *
# *                  2D Q4 Finite Element Analysis                             *
# *                     JIARUI WANG, ZIHAN WANG                                *
# *                         Copyright (c) 2026                                 *
# *                                                                            *
# ******************************************************************************
from fem2d_exact_solutions import PlateWithHole, CantileverBeam

# Boundary flags per (u_x, u_y): 0 = free, 1 = prescribed value, 99 = analytical.

problem_settings = {

    # Problem 1: plate with a hole
    'Plate_with_hole': {
        'solution': PlateWithHole,
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
        'dirichlet_values': {'FIXED_LEFT': (0.0, 0.0), 'FIXED_BOTTOM': (0.0, 0.0)},
        'traction': {'TRACTION_RIGHT': (99, 99), 'TRACTION_TOP': (99, 99)},
        'title': 'Plate with Hole',
        '2d_problem_type': 'plane_stress',
        'convergence': True
    },


    # Problem 2: cantilever beam
    'Cantilever_beam': {
        'solution': CantileverBeam,
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
        'traction': {'TRACTION_RIGHT': (99, 99)},
        'title': 'Cantilever Beam',
        '2d_problem_type': 'plane_strain',
        'convergence': True
    },


    # Problem 3: custom problem
    'Beam_with_tip_load': {
        'model_characteristics': [10.0, 2.0, 5.0],  # [length L, height H, tip load P]
        'mesh_file': '2d_beam_with_tip_load/model_cantilever_level_2.k',
        'material_properties': [1.0e3, 0.30],
        'integration_order': 2,
        'dirichlet': {'FIXED_LEFT': (1, 1)},
        'dirichlet_values': {'FIXED_LEFT': (0.0, 0.0)},
        'traction': {'TRACTION_RIGHT': (0, 1)},
        'traction_values': {'TRACTION_RIGHT': (0.0, 2.5)},
        'body_force': (0.0, 0.0),
        'title': 'Beam_with_tip_load',
        '2d_problem_type': 'plane_stress',
        'convergence': False
    },

}