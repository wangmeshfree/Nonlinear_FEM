# ******************************************************************************
# *                                                                            *
# *        MAE 5036 - Advanced Computational Solid Mechanics II                *
# *                                                                            *
# *                  2D Q4 Finite Element Analysis                             *
# *                     JIARUI WANG, ZIHAN WANG                                *
# *                         Copyright (c) 2026                                 *
# *                                                                            *
# ******************************************************************************


"""Analytical solutions: each class holds its parameters and gives u, sigma, body force, and traction."""

import numpy as np


class PlateWithHole:
    """Kirsch solution for an infinite plate with a traction-free circular hole.

    Quarter model x >= 0, y >= 0 under remote x-direction tension sigma0.
    theta is measured from +x; the material occupies r >= a. Rigid translation
    and rotation are zero, matching the quarter-model symmetry constraints.
    """

    def __init__(self, hole_radius, tension, E, nu, plane_problem_type):
        self.a = hole_radius
        self.sigma0 = tension
        self.E = E
        self.nu = nu
        self.plane_problem_type = plane_problem_type

    def displacement(self, x, y):
        a, sigma0, E, nu = self.a, self.sigma0, self.E, self.nu
        if self.plane_problem_type == 'plane_stress':
            beta = (3.0 - nu) / (1.0 + nu)
        elif self.plane_problem_type == 'plane_strain':
            beta = 3.0 - 4.0 * nu
        else:
            raise ValueError(f'Unknown plane problem type: {self.plane_problem_type}')

        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        mu = E / (2.0 * (1.0 + nu))
        ar = a / r
        scale = sigma0 * a / (8.0 * mu)
        ux = scale * (
            (r / a) * (beta + 1.0) * np.cos(theta)
            + 2.0 * ar * ((beta + 1.0) * np.cos(theta) + np.cos(3.0 * theta))
            - 2.0 * ar**3 * np.cos(3.0 * theta)
        )
        uy = scale * (
            (r / a) * (beta - 3.0) * np.sin(theta)
            + 2.0 * ar * ((1.0 - beta) * np.sin(theta) + np.sin(3.0 * theta))
            - 2.0 * ar**3 * np.sin(3.0 * theta)
        )
        return np.array([ux, uy])

    def stress(self, x, y):
        """Return [[sxx, sxy], [sxy, syy]]; points inside the hole are clamped to r = a."""
        a, sigma0 = self.a, self.sigma0
        r = max(np.sqrt(x**2 + y**2), a)
        th = np.arctan2(y, x)
        rr = a / r

        srr = 0.5*sigma0*(1.0 - rr**2.0) + 0.5*sigma0*(1.0 - 4.0*rr**2.0 + 3.0*rr**4.0)*np.cos(2.0*th)
        stt = 0.5*sigma0*(1.0 + rr**2.0) - 0.5*sigma0*(1.0 + 3.0*rr**4.0)*np.cos(2.0*th)
        srt = -0.5*sigma0*(1.0 + 2.0*rr**2.0 - 3.0*rr**4.0)*np.sin(2.0*th)

        RM = np.array([[np.cos(th), np.sin(th)],
                       [-np.sin(th), np.cos(th)]])
        return RM.T @ np.array([[srr, srt], [srt, stt]]) @ RM

    def body_force(self, x, y):
        return np.array([0.0, 0.0])

    def traction(self, x, y, normal):
        return self.stress(x, y) @ normal


class CantileverBeam:
    """Classical cantilever beam with unit thickness and a tip load.

    The beam occupies 0 <= x <= L and -H/2 <= y <= H/2. A positive tip load P
    acts in +y at x = L, and the Euler-Bernoulli stress field is
        sigma_xx = -M(x) * y / I,  sigma_xy = 3*P/(2*H) * (1 - 4*y^2/H^2),  sigma_yy = 0,
    with M(x) = P*(L - x) and I = H^3/12.
    """

    def __init__(self, length, height, tip_load, E, nu, plane_problem_type):
        self.length = length
        self.height = height
        self.tip_load = tip_load
        self.E = E
        self.nu = nu
        self.plane_problem_type = plane_problem_type

    def displacement(self, x, y):
        L, H, P = self.length, self.height, self.tip_load
        if self.plane_problem_type == 'plane_stress':
            E_eff, nu_eff = self.E, self.nu
        elif self.plane_problem_type == 'plane_strain':
            E_eff = self.E / (1.0 - self.nu**2)
            nu_eff = self.nu / (1.0 - self.nu)
        else:
            raise ValueError(f'Unknown plane problem type: {self.plane_problem_type}')

        second_moment = H**3 / 12.0
        ux = -P * y / (6.0 * E_eff * second_moment) * (
            (6.0 * L - 3.0 * x) * x
            + (2.0 + nu_eff) * (y**2 - H**2 / 4.0)
        )
        uy = P / (6.0 * E_eff * second_moment) * (
            3.0 * nu_eff * y**2 * (L - x)
            + (4.0 + 5.0 * nu_eff) * H**2 * x / 4.0
            + (3.0 * L - x) * x**2
        )
        return np.array([ux, uy])

    def stress(self, x, y):
        """Return [[sxx, sxy], [sxy, syy]]."""
        L, H, P = self.length, self.height, self.tip_load
        second_moment = H**3 / 12.0
        bending_moment = P * (L - x)

        sxx = -bending_moment * y / second_moment
        syy = 0.0
        sxy = (3.0 * P / (2.0 * H)) * (1.0 - 4.0 * y**2 / H**2)
        return np.array([[sxx, sxy], [sxy, syy]])

    def body_force(self, x, y):
        return np.array([0.0, 0.0])

    def traction(self, x, y, normal):
        return self.stress(x, y) @ normal
