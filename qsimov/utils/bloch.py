'''
QSimov: A Quantum Computing Toolkit.
Copyright (C) 2017  Hernán Indíbil de la Cruz Calvo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import matplotlib.pyplot as plt
import numpy as np


def plot_arc3d(vector1, vector2, radius=0.2, ax=None, colour="C0", color=None):
    """ Plot arc between two given vectors in 3D space. """
    '''
    https://stackoverflow.com/questions/47321839/how-to-show-the-angle-by-an-arc-between-two-3d-vectors-in-matplotlib
    answer by Phil G (https://stackoverflow.com/users/3912166/phil-g)
    Calculate vector between two vector end points, and the resulting
    spherical angles for various points along this vector. From this,
    derive points that lie along the arc between vector1 and vector2.
    '''
    if color is not None:
        colour = color
    v = [i - j for i, j in zip(vector1, vector2)]
    v_points_direct = [(vector2[0] + v[0] * num,
                        vector2[1] + v[1] * num,
                        vector2[2] + v[2] * num) for num in np.linspace(0, 1)]

    if any([np.linalg.norm(v_point) == 0.0 for v_point in v_points_direct]):
        return
    v_points_aux = [(v_point, np.linalg.norm(v_point))
                    for v_point in v_points_direct]
    v_thetas = [np.arccos(v_point[2]/norm)
                for v_point, norm in v_points_aux if norm != 0.0]
    v_phis = [np.arctan2(v_point[1], v_point[0])
              for v_point, norm in v_points_aux if norm != 0.0]

    v_points_arc = [(radius*np.sin(theta)*np.cos(phi),
                     radius*np.sin(theta)*np.sin(phi),
                     radius*np.cos(theta))
                    for theta, phi in zip(v_thetas, v_phis)]
    v_points_arc.append((0, 0, 0))

    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    """ Plot polygon.
    Face colour must be set afterwards, otherwise it over-rides
    the transparency.
    https://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
    """
    points_collection = Poly3DCollection([v_points_arc], alpha=0.4)
    points_collection.set_facecolor(colour)
    ax.add_collection3d(points_collection)


def draw_bloch_sphere(theta, phi, length=1.0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plt.axis("off")
    # ax.set_aspect("equal")  # NotImplementedError
    # Alternative to equal aspect:
    ax.set_box_aspect([ub - lb
                       for lb, ub in (getattr(ax, f"get_{a}lim")()
                                      for a in "xyz")])

    # draw axes
    # X axis
    ax.quiver(-1.1, 0, 0, 2.2, 0, 0, color="r",
              arrow_length_ratio=0.07, linewidths=0.75)
    # Y axis
    ax.quiver(0, -1.1, 0, 0, 2.2, 0, color="g",
              arrow_length_ratio=0.07, linewidths=0.75)
    # Z axis
    ax.quiver(0, 0, -1.1, 0, 0, 2.2, color="b",
              arrow_length_ratio=0.07, linewidths=0.75)

    # ax.plot([-1.1, 1.1],[0, 0], [0, 0],lw=1,color="r")  # X axis
    # ax.plot([0, 0],[-1.1, 1.1], [0, 0],lw=1,color="g")  # Y axis
    # ax.plot([0, 0], [0, 0],[-1.1, 1.1],lw=1,color="b")  # Z axis
    ax.text(-1.25, 0, 0, "|-\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.text(+1.25, 0, 0, "|+\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.text(0, -1.25, 0, "|-i\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.text(0, +1.25, 0, "|i\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.text(0, 0, -1.25, "|1\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.text(0, 0, +1.25, "|0\u27E9", fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    ax.scatter([0], [0], [0], color="k", s=10)  # Origin

    # draw sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color="white",
                    linewidth=0, alpha=0.25)

    # Calculate the qubit point
    x = np.cos(phi) * np.sin(theta) * length
    y = np.sin(phi) * np.sin(theta) * length
    z = np.cos(theta) * length
    # draw a point
    ax.scatter([x], [y], [z], color="k", s=10)
    '''
    coord_string = f"({round(x, 2)}, {round(y, 2)}, {round(z, 2)}"
    if length != 1.0:
        coord_string += f", {round(length, 2)})"
    coord_string += ")"
    ax.text(x * 1.25, y * 1.25, z * 1.25, coord_string, fontsize=12,
            horizontalalignment="center", verticalalignment="center")
    '''
    # draw a vector
    ax.quiver(0, 0, 0, x, y, z, color="k", linewidths=2)
    # draw dashed lines
    ax.plot([x, x], [y, y], [z, 0], lw=1, color="k", linestyle="dashed")
    ax.plot([x, 0], [y, 0], [0, 0], lw=1, color="k", linestyle="dashed")
    # draw angles
    # theta
    plot_arc3d([0, 0, 1], [x, y, z], radius=0.25*length, ax=ax, color="b")

    # phi
    if phi < np.pi:
        plot_arc3d([1, 0, 0], [x, y, 0], radius=0.25*length, ax=ax, color="r")
    elif phi >= np.pi:
        plot_arc3d([1, 0, 0], [0, 1, 0], radius=0.25*length, ax=ax, color="r")
        plot_arc3d([0, 1, 0], [-1, 0, 0], radius=0.25*length, ax=ax, color="r")
        plot_arc3d([-1, 0, 0], [x, y, 0], radius=0.25*length, ax=ax, color="r")

    ax.view_init(elev=15, azim=15)
    # plt.show()
    return fig


# draw_bloch_sphere(3*np.pi/4, 3.5*np.pi/2, length=1)
