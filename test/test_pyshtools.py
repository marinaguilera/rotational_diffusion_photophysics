from pyshtools import *
import spherical as spherical
import numpy as np
import matplotlib.pyplot as plt
from sympy import Matrix
from sympy.physics.quantum.cg import CG

lmax = 20
cos2 = SHGrid.from_zeros(lmax=lmax)
theta = cos2.lats() - 90
phi = cos2.lons()
omega = np.meshgrid(theta,phi)
omega[0] = omega[0]*np.pi/180
omega[1] = omega[1]*np.pi/180
probability = np.cos(omega[0])**2
cos2.data = probability.transpose()
cilm = cos2.expand()
vecclim = shio.SHCilmToVector(cilm.coeffs, lmax)

fig, ax = cos2.plot3d(cmap="Reds", scale=1e6, elevation=30, azimuth=90)
# ax.set_xlim3d(-1, 1)
# ax.set_ylim3d(-1, 1)
ax.set_zlim3d(-0.75, 0.75)

plt.figure()
plt.imshow(cilm.coeffs[0])
plt.figure()
plt.imshow(cilm.coeffs[1])



plt.show()


utils.Wigner3j(1,2,2,1,1)

cg = spherical.clebsch_gordan(1, 0, 0, 0, 3, 0)

M = Matrix([[CG(1,0,2,0,3,0), CG(1,0,3,0,3,0) ],[  CG(1,0,2,1,3,1), CG(1,0,2,2,2,2) ]])

M = Matrix([[3, -2,  4, -2],
            [5,  3, -3, -2],
            [5, -2,  2, -2],
            [5, -2, -3,  3]])
P, D = M.diagonalize()
