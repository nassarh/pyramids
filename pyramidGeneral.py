# Import stuff
import numpy as np
import pyramids as pd
import pyvista as pv  # for 3D rendering

# Define reference pyramid
# square basis, equilateral
u0 = np.array([1., 0., 0.])
v0 = np.array([0., 1., 0.])
wb0 = np.array([0.5, 0.5, 1 / np.sqrt(2)])
w0 = np.array([0.5, 0.5, -1 / np.sqrt(2)])

# opening angle
theta = 89 * np.pi / 180

# initial angle of incidence relative to axis of symmetry
beta = np.pi/2

# number of cells per parallel
N = 100

# number of cells per meridian
cells = 200

# Define zigzag: two vectors
u, v, _, _ = pd.zigcircle(theta, beta, N, u0=u0, v0=v0, w0=w0, wb0=wb0)

# define rotational symmetry by a general rotation
# axis
p = np.array([1, 1, 1])
p = p / np.linalg.norm(p)
P = np.array([[0, -p[2], p[1]], [p[2], 0, -p[0]], [-p[1], p[0], 0]])

# angle
dphi = 2 * np.pi / N

# rotation
rot = np.eye(3) + np.sin(dphi) * P + (1 - np.cos(dphi)) * np.dot(P, P)

# Build pattern: crease vectors and vertices
U, V, W, _ = pd.manysteps(u, v, u0, v0, w0, wb0, cells, rot=rot)
X, Y, Z, Xt, Yt, Zt = pd.integrate(U, V, W, N, per=False, rot=rot)

x = np.hstack([X.ravel(),Xt.ravel()])
y = np.hstack([Y.ravel(),Yt.ravel()])
z = np.hstack([Z.ravel(),Zt.ravel()])

nodes = np.hstack([x[:,None],y[:,None],z[:,None]])

# plot
pl = pv.Plotter()
pl.set_background("black", top="white")
grid = pv.make_tri_mesh(nodes,np.array(pd.triangles(N,cells)))
pl.add_mesh(grid, show_edges=True, line_width=1, color=[193, 225, 193])
pl.show_axes()
pl.show()