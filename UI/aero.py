"""Wall and field post-processing of CMFSN structured-grid solutions (.vts read with pyvista).

Shared by visualize.py and grid_convergence.py. Every function takes a pyvista StructuredGrid
whose j = 0 row of cells touches the wall (O-mesh) and whose cell data hold the exported fields.
"""
from __future__ import annotations

import numpy as np


def _field(mesh, name, default):
    return float(mesh.field_data[name][0]) if name in mesh.field_data else default


WALL_PRESSURE_METHODS = ("cell", "extrapolated")
# "cell" is the default: it is the pressure of the solver's own wall flux. The slip wall copies the first
# cell's pressure into the ghost cell (zero normal pressure gradient), so the face pressure is p1. On a
# curved wall the true gradient is rho V^2 / R; the solver's first rows do not follow it, and a linear
# extrapolation from them moves the coefficients away from the discrete solution instead of improving them.


def wall_pressure(mesh, method="cell"):
    """Pressure on every wall face (j = 0), in the order of the wall cells.

    "cell": the pressure of the first cell, at its centre (first order in the wall distance).
    "extrapolated": linear extrapolation to the face along the wall normal from the first two cells,
    p_w = p1 + (p1 - p2) * d1 / (d2 - d1), with d1, d2 their centres' distances to the face line.
    Falls back to "cell" where the mesh has a single row of cells or the two distances are not increasing.
    """
    if method not in WALL_PRESSURE_METHODS:
        raise ValueError(f"unknown wall pressure method {method!r}, expected one of {WALL_PRESSURE_METHODS}")
    ni, nj, _ = mesh.dimensions
    pressure = np.asarray(mesh.cell_data["pressure"], dtype=float).reshape(nj - 1, ni - 1)
    if method == "cell" or nj < 3:
        return pressure[0].copy()
    nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:3, :, :2]
    middle = .5 * (nodes[0, :-1] + nodes[0, 1:])
    edge = nodes[0, 1:] - nodes[0, :-1]
    unit = np.column_stack([edge[:, 1], -edge[:, 0]]) / np.linalg.norm(edge, axis=1)[:, None]
    centre = lambda j: .25 * (nodes[j, :-1] + nodes[j, 1:] + nodes[j + 1, :-1] + nodes[j + 1, 1:])
    d1 = np.abs(np.einsum("ij,ij->i", centre(0) - middle, unit))
    d2 = np.abs(np.einsum("ij,ij->i", centre(1) - middle, unit))
    valid = d2 - d1 > 1e-12 * np.maximum(d2, 1e-300)
    result = pressure[0].copy()
    result[valid] += (pressure[0, valid] - pressure[1, valid]) * d1[valid] / (d2[valid] - d1[valid])
    return result


def aerodynamic_coefficients(mesh, chord_ref=None, x_ref=None, pressure_method="cell"):
    """Force coefficients on the wall (j = 0).

    Pressure: wall_pressure(mesh, pressure_method), the first cell's by default. The friction
    part is present only when the solution stores a Reynolds number: tau = mu * tangential velocity /
    wall distance, from the first cell centres, with the exported Sutherland viscosity of the first
    cells (the wall ghost cells share their temperature at an adiabatic wall), or the constant
    mu_inf = M_inf * sqrt(gamma) / Re for files written before the viscosity was exported.

    Reference length: chord_ref, or the measured chord (x extent of the wall) when None.
    Moment centre: (x_ref, y of the leading edge), or the quarter of the measured chord when None.
    Cm is positive nose-up. Returns None for a free stream at rest.
    """
    mach, gamma, alpha = _field(mesh, "mach_inf", 0.0), _field(mesh, "gamma", 1.4), _field(mesh, "alpha_rad", 0.0)
    rho_inf, p_inf, reynolds = _field(mesh, "rho_inf", 1.0), _field(mesh, "pressure_inf", 1.0), _field(mesh, "reynolds_inf", np.inf)
    dynamic_pressure = .5 * gamma * p_inf * mach * mach
    ni, nj, _ = mesh.dimensions
    if not dynamic_pressure > 0 or ni < 2 or nj < 2:
        return None
    points = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:2, :, :2]
    middle = .5 * (points[0, :-1] + points[0, 1:])
    centers = .25 * (points[0, :-1] + points[0, 1:] + points[1, :-1] + points[1, 1:])
    edge = points[0, 1:] - points[0, :-1]
    length = np.linalg.norm(edge, axis=1)
    normal = np.column_stack([edge[:, 1], -edge[:, 0]])  # Scaled by the face length.
    normal *= np.sign(np.einsum("ij,ij->i", normal, middle - centers))[:, None]  # From the fluid into the body.
    row = lambda name: np.asarray(mesh.cell_data[name], dtype=float).reshape(nj - 1, ni - 1, -1)[0]
    pressure_force = (wall_pressure(mesh, pressure_method) - p_inf)[:, None] * normal
    friction_force = np.zeros_like(pressure_force)
    viscous = bool(np.isfinite(reynolds) and reynolds > 0)
    if viscous:
        distance = np.abs(np.einsum("ij,ij->i", normal, middle - centers))  # Wall distance times face length.
        valid = distance > 0
        unit = normal[valid] / length[valid, None]
        velocity = row("velocity")[valid, :2]
        tangential = velocity - np.einsum("ij,ij->i", velocity, unit)[:, None] * unit
        if "viscosity" in mesh.cell_data:
            viscosity = row("viscosity")[valid, 0]
        else:
            viscosity = np.full(int(valid.sum()), rho_inf * mach * np.sqrt(gamma * p_inf / rho_inf) / reynolds)
        friction_force[valid] = viscosity[:, None] * tangential * (length[valid] ** 2 / distance[valid])[:, None]
    wall_x = points[0, :, 0]
    measured_chord = float(np.ptp(wall_x))
    chord = measured_chord if chord_ref is None else float(chord_ref)
    if not chord > 0:
        return None
    leading_edge = points[0, np.argmin(wall_x)]
    origin = np.array([leading_edge[0] + .25 * measured_chord if x_ref is None else float(x_ref), leading_edge[1]])
    lever = middle - origin

    def coefficients(force):
        fx, fy = np.nansum(force, axis=0) / (dynamic_pressure * chord)
        moment = -np.nansum(lever[:, 0] * force[:, 1] - lever[:, 1] * force[:, 0]) / (dynamic_pressure * chord * chord)
        return np.array([-fx * np.sin(alpha) + fy * np.cos(alpha), fx * np.cos(alpha) + fy * np.sin(alpha), moment])

    pressure_part, friction_part = coefficients(pressure_force), coefficients(friction_force)
    total = pressure_part + friction_part
    return {"cl": total[0], "cd": total[1], "cm": total[2], "cd_pressure": pressure_part[1], "cd_friction": friction_part[1], "viscous": viscous, "chord": chord}


def wall_distribution(mesh, pressure_method="cell"):
    """Wall Cp against x (pressure when the free stream is at rest), from wall_pressure(mesh, pressure_method).

    With "extrapolated", x is the middle of each wall face; with "cell", the first cell centre.
    """
    mach, gamma, p_inf = _field(mesh, "mach_inf", 0.0), _field(mesh, "gamma", 1.4), _field(mesh, "pressure_inf", 1.0)
    dynamic_pressure = .5 * gamma * p_inf * mach * mach
    ni, nj, _ = mesh.dimensions
    nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:2, :, 0]
    if pressure_method == "cell":
        x = .25 * (nodes[0, :-1] + nodes[0, 1:] + nodes[1, :-1] + nodes[1, 1:])
    else:
        x = .5 * (nodes[0, :-1] + nodes[0, 1:])
    pressure = wall_pressure(mesh, pressure_method)
    if dynamic_pressure > 0:
        return "cp", x, (pressure - p_inf) / dynamic_pressure
    return "pressure", x, pressure


def cell_areas(mesh):
    """Area of every quadrilateral cell (shoelace on its four nodes), ordered like the cell data."""
    ni, nj, _ = mesh.dimensions
    p = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:, :, :2]
    a, b, c, d = p[:-1, :-1], p[:-1, 1:], p[1:, 1:], p[1:, :-1]
    cross = lambda u, v: u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    return np.abs(.5 * (cross(a, b) + cross(b, c) + cross(c, d) + cross(d, a))).ravel()


def entropy_error(mesh, exclude_radius=0.0):
    """Area-weighted L2 norm of (s - s_inf) / s_inf over the cells, with s = p / rho^gamma.

    In a smooth inviscid (Euler) flow without shocks the entropy is uniform, so this measures the
    discretization error; its slope against the cell size gives the spatial order of the scheme.
    exclude_radius > 0 leaves out the cells whose centre lies within that distance of the trailing
    edge (wall node of largest x), a singular point that can lower the measured order.
    """
    gamma, rho_inf, p_inf = _field(mesh, "gamma", 1.4), _field(mesh, "rho_inf", 1.0), _field(mesh, "pressure_inf", 1.0)
    s_inf = p_inf / rho_inf ** gamma
    pressure = np.asarray(mesh.cell_data["pressure"], dtype=float).ravel()
    density = np.asarray(mesh.cell_data["density"], dtype=float).ravel()
    error = (pressure / density ** gamma - s_inf) / s_inf
    areas = cell_areas(mesh)
    keep = np.ones(error.size, dtype=bool)
    if exclude_radius > 0:
        ni, nj, _ = mesh.dimensions
        nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:, :, :2]
        trailing_edge = nodes[0, np.argmax(nodes[0, :, 0])]
        centres = (.25 * (nodes[:-1, :-1] + nodes[:-1, 1:] + nodes[1:, 1:] + nodes[1:, :-1])).reshape(-1, 2)
        keep = np.linalg.norm(centres - trailing_edge, axis=1) > exclude_radius
    return float(np.sqrt(np.sum(areas[keep] * error[keep] ** 2) / np.sum(areas[keep])))
