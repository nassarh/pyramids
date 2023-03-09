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
theta = 90 * np.pi / 180

# number of cells per parallel
N = 5

# number of cells per meridian
cells = 5

# Define zigzag: two vectors + invariance by translation (default)
u, v, _ = pd.zigzag(theta, u0=u0, v0=v0, w0=w0, wb0=wb0)

# Build pattern: crease vectors and vertices
U, V, W, _ = pd.manysteps(u, v, u0, v0, w0, wb0, cells)
X, Y, Z, Xt, Yt, Zt = pd.integrate(U, V, W, N, per=False)

x = np.hstack([X.ravel(),Xt.ravel()])
y = np.hstack([Y.ravel(),Yt.ravel()])
z = np.hstack([Z.ravel(),Zt.ravel()])

nodes = np.hstack([x[:,None],y[:,None],z[:,None]])


# plot
pl = pv.Plotter()
pl.set_background("black", top="white")
grid = pv.make_tri_mesh(nodes,np.array(pd.triangles(N,cells)))
pl.add_mesh(grid, show_edges=True, line_width=3, color=[193, 225, 193])
pl.camera_position = 'xy'
pl.show_axes()
pl.show()
