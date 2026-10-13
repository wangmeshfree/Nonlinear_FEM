"""Three-dimensional analytical benchmarks (no plane-stress/strain reduction)."""
import numpy as np


def elasticity(E, nu):
    """3D Hooke matrix for [xx, yy, zz, yz, xz, xy] engineering strains."""
    if E <= 0.0:
        raise ValueError('Young modulus must be positive')
    if nu <= -1.0 or nu >= 0.5:
        raise ValueError('Poisson ratio must satisfy -1 < nu < 0.5')

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


def stress_tensor(s):
    # Voigt order: xx, yy, zz, yz, xz, xy; engineering shear strains.
    return np.array([[s[0],s[5],s[4]],[s[5],s[1],s[3]],[s[4],s[3],s[2]]])


class UniformTension:
    """Exact affine displacement under remote x tension in a rectangular solid."""
    def __init__(self, tension, E, nu):
        self.tension,self.E,self.nu=tension,E,nu
    def displacement(self,x,y,z):
        return self.tension/self.E*np.array([x,-self.nu*y,-self.nu*z])
    def stress(self,x,y,z):
        return np.diag([self.tension,0.,0.])
    def body_force(self,x,y,z):return np.zeros(3)
    def traction(self,x,y,z,normal):return self.stress(x,y,z)@normal


class ManufacturedSine:
    """Exact u=(sin(pi*x)sin(pi*y)sin(pi*z),0,0) on the unit cube.

    Body force is derived from b=-div(sigma); all six faces have zero displacement.
    """
    def __init__(self, amplitude, E, nu):
        self.A=amplitude;self.E=E;self.nu=nu;self.D=elasticity(E,nu)
        self.mu=E/(2*(1+nu));self.lam=E*nu/((1+nu)*(1-2*nu))
    def displacement(self,x,y,z):
        return np.array([self.A*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z),0.,0.])
    def stress(self,x,y,z):
        sx,sy,sz=np.sin(np.pi*np.array([x,y,z]))
        cx,cy,cz=np.cos(np.pi*np.array([x,y,z]))
        strain=self.A*np.pi*np.array([cx*sy*sz,0.,0.,0.,sx*sy*cz,sx*cy*sz])
        return stress_tensor(self.D@strain)
    def body_force(self,x,y,z):
        sx,sy,sz=np.sin(np.pi*np.array([x,y,z]))
        cx,cy,cz=np.cos(np.pi*np.array([x,y,z]))
        return self.A*np.pi**2*np.array([(self.lam+4*self.mu)*sx*sy*sz,
                         -(self.lam+self.mu)*cx*cy*sz,-(self.lam+self.mu)*cx*sy*cz])
    def traction(self,x,y,z,normal):return self.stress(x,y,z)@normal


class ManufacturedZeroTraction:
    """Unit cube: prescribed displacement on x faces, zero traction on y/z faces.

    u_x = A * 4*x*(1-x) * 16*y^2*(1-y)^2 * 16*z^2*(1-z)^2.
    u_y = u_z = 0. Both displacement and its gradient vanish on y/z faces.
    """
    def __init__(self, amplitude, E, nu):
        self.A = amplitude
        self.E = E
        self.nu = nu
        self.mu = E / (2.0 * (1.0 + nu))
        self.lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        self.D = elasticity(E, nu)

    def displacement(self, x, y, z):
        fx = 4.0 * x * (1.0 - x)
        fy = 16.0 * y**2 * (1.0 - y)**2
        fz = 16.0 * z**2 * (1.0 - z)**2
        ux = self.A * fx * fy * fz
        return np.array([ux, 0.0, 0.0])

    def stress(self, x, y, z):
        fx = 4.0 * x * (1.0 - x)
        fy = 16.0 * y**2 * (1.0 - y)**2
        fz = 16.0 * z**2 * (1.0 - z)**2
        dfx = 4.0 - 8.0 * x
        dfy = 32.0 * y * (1.0 - y) * (1.0 - 2.0 * y)
        dfz = 32.0 * z * (1.0 - z) * (1.0 - 2.0 * z)

        dux_dx = self.A * dfx * fy * fz
        dux_dy = self.A * fx * dfy * fz
        dux_dz = self.A * fx * fy * dfz
        strain = np.array([dux_dx, 0.0, 0.0, 0.0, dux_dz, dux_dy])
        stress_vector = self.D @ strain
        return stress_tensor(stress_vector)

    def body_force(self, x, y, z):
        fx = 4.0 * x * (1.0 - x)
        fy = 16.0 * y**2 * (1.0 - y)**2
        fz = 16.0 * z**2 * (1.0 - z)**2
        dfx = 4.0 - 8.0 * x
        dfy = 32.0 * y * (1.0 - y) * (1.0 - 2.0 * y)
        dfz = 32.0 * z * (1.0 - z) * (1.0 - 2.0 * z)
        ddfx = -8.0
        ddfy = 32.0 * (1.0 - 6.0 * y + 6.0 * y**2)
        ddfz = 32.0 * (1.0 - 6.0 * z + 6.0 * z**2)

        dux_dxx = self.A * ddfx * fy * fz
        dux_dyy = self.A * fx * ddfy * fz
        dux_dzz = self.A * fx * fy * ddfz
        dux_dxy = self.A * dfx * dfy * fz
        dux_dxz = self.A * dfx * fy * dfz

        bx = -(self.lam + 2.0 * self.mu) * dux_dxx
        bx -= self.mu * dux_dyy
        bx -= self.mu * dux_dzz
        by = -(self.lam + self.mu) * dux_dxy
        bz = -(self.lam + self.mu) * dux_dxz
        return np.array([bx, by, bz])

    def traction(self, x, y, z, normal):
        stress = self.stress(x, y, z)
        traction = stress @ normal
        return traction
