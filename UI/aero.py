"""Wall and field post-processing of CMFSN structured-grid solutions (.vts read with pyvista).

Shared by visualize.py and grid_convergence.py. Every function takes a pyvista StructuredGrid
whose j = 0 row of cells touches the wall (O-mesh) and whose cell data hold the exported fields.
"""
from __future__ import annotations

import numpy as np


def _field(mesh, name, default):
    return float(mesh.field_data[name][0]) if name in mesh.field_data else default


WALL_PRESSURE_METHODS = ("auto", "rizzi", "cell", "extrapolated")
# "auto" is the default: the pressure of the solver's own wall flux.
#   Euler2D (no "reynolds_inf" in the file): the slip wall extrapolates the pressure into the ghost cells with
#   Rizzi's normal-momentum relation and imposes P_f = p1 - P_eta / 2 on the wall face -> "rizzi".
#   RANS2D: the no-slip ghost copies p1 (zero normal pressure gradient), so the face pressure is p1 -> "cell".
# "extrapolated" (linear, from the first two cells) is kept for comparison only.


def _is_viscous(mesh):
    reynolds = _field(mesh, "reynolds_inf", np.inf)
    return bool(np.isfinite(reynolds) and reynolds > 0)


def resolve_wall_pressure_method(mesh, method):
    """The method actually used: "auto" becomes "rizzi" for an inviscid solution, "cell" for a viscous one."""
    if method not in WALL_PRESSURE_METHODS:
        raise ValueError(f"unknown wall pressure method {method!r}, expected one of {WALL_PRESSURE_METHODS}")
    if method == "auto":
        return "cell" if _is_viscous(mesh) else "rizzi"
    return method


def _cell_centroids(nodes):
    """Centroid of every quadrilateral cell, nodes of shape (nj, ni, 2) -> (nj - 1, ni - 1, 2).
    Same point as the solver's eval_CellCenterAreas (area-weighted centroid of the cell)."""
    corners = (nodes[:-1, :-1], nodes[:-1, 1:], nodes[1:, 1:], nodes[1:, :-1])
    twice_area = np.zeros(corners[0].shape[:-1])
    moment = np.zeros(corners[0].shape)
    for a in range(4):
        p, q = corners[a], corners[(a + 1) % 4]
        cross = p[..., 0] * q[..., 1] - q[..., 0] * p[..., 1]
        twice_area += cross
        moment += (p + q) * cross[..., None]
    return moment / (3 * twice_area[..., None])


def _rizzi_wall_pressure(mesh):
    """Face pressure of the Euler2D slip wall, P_f = p1 - P_eta / 2, with the solver's formula (CellOps.h,
    ApplyWallConditions and eval_wall_metrics):
        P_eta = [(x_eps . x_eta) p_eps + rho1 det(V1, x_eta) det(x_eps_eps, V1)] / |x_eps|^2
    x_eps, x_eps_eps: central differences of the wall face midpoints; x_eta: from the first ghost centre
    (mirror of the wall cell centre through the face midpoint) to the wall cell centre; p_eps: central
    difference of the wall cell pressures. Periodic in i, like the solver: closed O-mesh."""
    ni, nj, _ = mesh.dimensions
    nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:2, :, :2]
    first_row = lambda name: np.asarray(mesh.cell_data[name], dtype=float).reshape(nj - 1, ni - 1, -1)[0]
    p1, rho1, v1 = first_row("pressure")[:, 0], first_row("density")[:, 0], first_row("velocity")[:, :2]
    middle = .5 * (nodes[0, :-1] + nodes[0, 1:])
    right, left = np.roll(middle, -1, axis=0), np.roll(middle, 1, axis=0)
    x_eps = .5 * (right - left)
    x_eps_eps = right - 2 * middle + left
    x_eta = 2 * (_cell_centroids(nodes)[0] - middle)
    p_eps = .5 * (np.roll(p1, -1) - np.roll(p1, 1))
    det = lambda u, v: u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    p_eta = (np.einsum("ij,ij->i", x_eps, x_eta) * p_eps + rho1 * det(v1, x_eta) * det(x_eps_eps, v1)) / np.einsum("ij,ij->i", x_eps, x_eps)
    return p1 - .5 * p_eta


def wall_pressure(mesh, method="auto"):
    """Pressure on every wall face (j = 0), in the order of the wall cells.

    "auto": "rizzi" for an inviscid solution, "cell" for a viscous one (the solver's wall flux).
    "rizzi": the Euler2D wall face pressure P_f = p1 - P_eta / 2 (see _rizzi_wall_pressure).
    "cell": the pressure of the first cell, at its centre.
    "extrapolated": linear extrapolation to the face along the wall normal from the first two cells,
    p_w = p1 + (p1 - p2) * d1 / (d2 - d1), with d1, d2 their centres' distances to the face line.
    Falls back to "cell" where the mesh has a single row of cells or the two distances are not increasing.
    """
    method = resolve_wall_pressure_method(mesh, method)
    ni, nj, _ = mesh.dimensions
    if method == "rizzi":
        return _rizzi_wall_pressure(mesh)
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


def aerodynamic_coefficients(mesh, chord_ref=None, x_ref=None, pressure_method="auto"):
    """Force coefficients on the wall (j = 0).

    Pressure: wall_pressure(mesh, pressure_method), the solver's wall face pressure by default. The friction
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
    viscous = _is_viscous(mesh)
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


def wall_distribution(mesh, pressure_method="auto"):
    """Wall Cp against x (pressure when the free stream is at rest), from wall_pressure(mesh, pressure_method).

    With "rizzi" or "extrapolated", x is the middle of each wall face; with "cell", the first cell centre.
    """
    mach, gamma, p_inf = _field(mesh, "mach_inf", 0.0), _field(mesh, "gamma", 1.4), _field(mesh, "pressure_inf", 1.0)
    dynamic_pressure = .5 * gamma * p_inf * mach * mach
    ni, nj, _ = mesh.dimensions
    nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:2, :, 0]
    pressure_method = resolve_wall_pressure_method(mesh, pressure_method)
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
