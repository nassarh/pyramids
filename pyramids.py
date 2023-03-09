# by Hussein Nassar (nassarh@missouri.edu)
# you are free to use and modify for research and education purposes with
# proper citation and attribution.
# for commercial use, bugs and questions, contact the author.

import numpy as np
import numpy.linalg as la


#########################################
# Tools to construct intersection of three spheres
#########################################
def w(x, y, x0, y0, z0, sgn=1):
    # Find intersection z of three spheres centered at 0, x and y
    # and of radius vectors z0, z0-x0 and z0-y0
    # (x,y,z) is right-handed if sgn > 0

    # Define basis (xs,ys,n) dual to (x,y,n)
    # n unitary orthogonal to (x,y)
    # (x,y,n) is right handed
    n = np.cross(x, y)
    J = la.norm(n)

    # Check if basis is degenerate
    if J == 0.0:
        raise ValueError("basis degenerate - angle too small")

    n = n / J
    xs = np.cross(y, n) / J
    ys = np.cross(n, x) / J

    # Define solution coordinates
    z1 = np.dot(z0, x0)
    z2 = np.dot(z0, y0)

    z3squared = (
        la.norm(z0) ** 2
        - la.norm(z1 * xs) ** 2
        - la.norm(z2 * ys) ** 2
        - 2 * z1 * z2 * np.dot(xs, ys)
    )

    # Check if solution exists
    if z3squared < 0:
        raise ValueError(
            "no solution - angle too large (well, it could be too small too)"
        )

    z3 = np.sqrt(z3squared)

    # Build intersection with given orientation (sgn)
    z = z1 * xs + z2 * ys + np.sign(sgn) * z3 * n

    return z


#########################################
# Cauchy's (initial value) problem under periodic boundary conditions
#########################################
def onestep(
    u,
    v,
    u0,
    v0,
    w0,
    wb0,
    rot=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
):
    # Apply the three spheres intersection for one dent twice to complete one pyramid
    # Takes into account periodic bc through rot
    # roll = 1 => it is a "principal" meridian
    #       -1 => it is a "secondary" meridian

    # Setup:
    #                                *                              *
    #                              * * *                          *
    #                           *    *    *                    *
    #                rot-1(UU)      Wb      VV              UU
    #                     *          *          *         *
    #                  *             *             *   *
    #  rot-1 <-      *********************************      -> rot
    #                  *             *             *
    #                     *          *          *
    #                        V      WW       U
    #                           *    *    *
    #                              * * *
    #                                *

    ww = w(u, v, u0, v0, w0)
    wwb = w(u - ww, v - ww, u0 - w0, v0 - w0, wb0)
    uu = np.dot(rot, wwb + ww - v)
    vv = wwb + ww - u

    return ww, wwb, uu, vv


def manysteps(
    u,
    v,
    u0,
    v0,
    w0,
    wb0,
    cells,
    rot=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
):
    # Iterate onestep for many cells
    # flags if construction was interrupted

    # Initialization
    U = np.zeros((cells + 1, 3))
    V = np.zeros((cells + 1, 3))
    W = np.zeros((2 * cells, 3))

    U[0, :] = u
    V[0, :] = v

    # Be optimistic
    flag = 0

    for n in range(0, cells):
        try:
            ww, wwb, uu, vv = onestep(U[n, :], V[n, :], u0, v0, w0, wb0, rot=rot)
            W[2 * n, :] = ww
            W[2 * n + 1, :] = wwb

            U[n + 1, :] = uu
            V[n + 1, :] = vv

        # Stop construction when the paper tears or self-contact
        except ValueError as e:
            print(e)
            print("after " + str(n) + " iterations")
            # make sure to keep a whole number of unit cells
            V = V[: n + 1, :]
            U = U[: n + 1, :]
            W = W[: 2 * n, :]
            # You can't always get what you want - MJ
            flag = 1
            break

    return U, V, W, flag


def integrate(
    U,
    V,
    W,
    N,
    rot=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    per=False,
):
    # Construct surface of revolution from fields U and V and rotation rot
    # N is the number of unit cells to be constructed along a parallel
    # per = False : first and last meridian are not identical
    # per = True  : first and last meridian are identical
    # retruns positions of basis nodes and pyramid tips

    # Initialization
    c = np.shape(U)[0]
    X = np.zeros((N + 1, c + 1))
    Y = np.zeros((N + 1, c + 1))
    Z = np.zeros((N + 1, c + 1))

    Xt = np.zeros((N, c - 1))
    Yt = np.zeros((N, c - 1))
    Zt = np.zeros((N, c - 1))

    # Integrate in the meridian direction
    for j in range(1, c + 1):
        [X[0, j], Y[0, j], Z[0, j]] = [X[0, j - 1], Y[0, j - 1], Z[0, j - 1]] + U[
            j - 1, :
        ]

    # Initialize current (dynamic) meridian
    UU = np.copy(U[0, :])
    VV = np.copy(V)
    WW = np.copy(W[::2, :])

    # Integrate in the parallel direction
    for i in range(1, N + 1):
        [X[i, 1:], Y[i, 1:], Z[i, 1:]] = [
            X[i - 1, :-1],
            Y[i - 1, :-1],
            Z[i - 1, :-1],
        ] + VV[:, :].T

        # Update U
        UU = np.dot(UU, rot)

        [X[i, 0], Y[i, 0], Z[i, 0]] = (
            [X[i - 1, 0], Y[i - 1, 0], Z[i - 1, 0]] + VV[0, :] - UU
        )

        # Update V
        VV = np.dot(VV, rot)

        # Find pyramid tips
        [Xt[i - 1, :], Yt[i - 1, :], Zt[i - 1, :]] = [
            X[i - 1, :-2],
            Y[i - 1, :-2],
            Z[i - 1, :-2],
        ] + WW.T

        # Update W
        WW = np.dot(WW, rot)

    if per:
        # Make sure ends meet
        [X[N, :], Y[N, :], Z[N, :]] = [X[0, :], Y[0, :], Z[0, :]]

    return X, Y, Z, Xt, Yt, Zt


#########################################
# Tools to define initial conditions
#########################################
def zigzag(
    theta,
    u0=np.array([1.0, 0.0, 0.0]),
    v0=np.array([-1.0, 0.0, 0.0]),
    w0=np.array([1.0, 0.0, 0.0]),
    wb0=np.array([1.0, 0.0, 0.0]),
):
    # Defines a uniform boundary condition
    # characterized by a single dent (u,v) with an opening theta
    # parallel direction is (1,0,0)
    # tangent plane is normal to (0,0,1)
    # flags if initial condition is not feasible

    # this is ok ...
    u = la.norm(u0) * np.array([np.sin(theta / 2), np.cos(theta / 2), 0])
    v = la.norm(v0) * np.array([-np.sin(theta / 2), np.cos(theta / 2), 0])

    # ... but let's adjust the parallel direction and the tangent plane

    # Initialize unit cell
    U, V, W, flag = manysteps(u, v, u0, v0, w0, wb0, 1)

    # Test if initial condition is feasible
    if flag:
        return u0, v0, flag

    # Define current, normalized, local basis
    t1 = u - v
    t1 = t1 / la.norm(t1)

    t2 = W[0, :] + W[1, :]
    t2 = t2 / la.norm(t2)

    n = np.cross(t1, t2)
    n = n / la.norm(n)

    # Complete the basis (t1, ?, n)
    t1p = np.cross(n, t1)

    # Adjust the dent
    u = np.array([np.dot(u, t1), np.dot(u, t1p), np.dot(u, n)])
    v = np.array([np.dot(v, t1), np.dot(v, t1p), np.dot(v, n)])

    return u, v, 0


def zigcircle(
    theta,
    beta,
    N=10,
    u0=np.array([1.0, 0.0, 0.0]),
    v0=np.array([-1.0, 0.0, 0.0]),
    w0=np.array([1.0, 0.0, 0.0]),
    wb0=np.array([1.0, 0.0, 0.0]),
):
    # Defines an axisymmetric boundary condition
    # characterized by a single dent (u,v) with an opening theta
    # N = nb of unit cells per parallel
    # beta = inclination of the meridian with respect to the axis of symmetry
    # returns the dent (u,v) and a rotation rot definining the symmetry
    # parallel is tangent to (1,0,0) at (0,0,0)
    # axis of symmetry is parallel to (0,1,0)
    # flags if initial condition is not feasible

    # Initialize dent
    u, v, flag = zigzag(theta, u0=u0, v0=v0, w0=w0, wb0=wb0)

    # Define local basis
    t1 = np.array([1.0, 0.0, 0.0])

    # Rotate dent by beta around t1
    u = (
        np.cos(beta) * (u - np.dot(u, t1) * t1)
        + np.sin(beta) * np.cross(t1, u)
        + np.dot(u, t1) * t1
    )
    v = (
        np.cos(beta) * (v - np.dot(v, t1) * t1)
        + np.sin(beta) * np.cross(t1, v)
        + np.dot(v, t1) * t1
    )

    # Define the symmetry of the pattern: rotation by 2pi/N about (0.,1.,0.)
    c = np.cos(2 * np.pi / N)
    s = np.sin(2 * np.pi / N)
    rot = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])

    return u, v, rot, flag


#########################################
# tools for plotting
#########################################
def triangles(N, c):
    # returns a list of triangles [p q r] that correspond to the facets
    # of a pyramid tessellation with N meridians and c pyramids per meridian
    # the numbering is the one dictated by pyramids.integrate for the basis nodes
    # followed by the tips

    nodes = (N + 1) * (c + 2)

    tr = (
        [
            [nodes - 2 * j + i + (c + 2) * j, i + (c + 2) * j, i + (c + 2) * j + 1]
            for i in range(c)
            for j in range(N)
        ]
        + [
            [
                nodes - 2 * j + i + (c + 2) * j,
                i + (c + 2) * j + 1,
                i + (c + 2) * j + c + 4,
            ]
            for i in range(c)
            for j in range(N)
        ]
        + [
            [
                nodes - 2 * j + i + (c + 2) * j,
                i + (c + 2) * j + c + 4,
                i + (c + 2) * j + c + 3,
            ]
            for i in range(c)
            for j in range(N)
        ]
        + [
            [nodes - 2 * j + i + (c + 2) * j, i + (c + 2) * j + c + 3, i + (c + 2) * j]
            for i in range(c)
            for j in range(N)
        ]
    )

    return tr


# Test stuff locally
if __name__ == "__main__":
    print(triangles(2, 5))
