"""Grid study post-processing: coefficients against Vassberg & Jameson (J. Aircraft 47(4), 2010),
residual histories, wall Cp, Mach lines, and the spatial order of the scheme from the entropy error.

Usage:   python UI/grid_convergence.py                       (most recent results/study_*)
         python UI/grid_convergence.py results/study_0001 --exclude-radius 0.05 --format pdf
Output:  <study>/post/ : table.csv, residuals, coefficients, entropy and one solution figure per mesh.

The study directory is written by Main (RunStudy): study.csv, then one sub-directory per mesh with
solution.vts and residuals.csv. Coefficients use the paper's reference chord 1 and moment centre
x = 0.25 unless --chord-ref / --x-ref say otherwise.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pyvista as pv  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))  # aero.py, next to this file.
from aero import WALL_PRESSURE_METHODS, aerodynamic_coefficients, entropy_error, resolve_wall_pressure_method, wall_distribution  # noqa: E402

WALL_PRESSURE_LABELS = {"rizzi": "Rizzi wall face pressure", "cell": "first cell centres", "extrapolated": "extrapolated to the wall"}


# Per-mesh coefficients of the paper, keyed by (Mach, alpha in degrees) then NC: (Cl, Cd, Cm).
# FLO82 (H-CUSP, point vortex at the far field): Tables 2-3. CFL3D without point vortex: Tables 6 and 8,
# closer to this solver, whose Riemann far field has no point vortex either. At alpha = 0, Cl = Cm = 0.
FLO82 = {
    (0.5, 0.0): {32: (0, .001979418, 0), 64: (0, .002270216, 0), 128: (0, .000307646, 0), 256: (0, -.000013299, 0),
                 512: (0, -.000015191, 0), 1024: (0, -.000005076, 0), 2048: (0, -.000001404, 0), 4096: (0, -.000000365, 0)},
    (0.5, 1.25): {32: (.181639135, .002482085, -.002926241), 64: (.178822451, .001974608, -.002487849),
                  128: (.180014345, .000298397, -.002294411), 256: (.180458183, -.000011045, -.002301404),
                  512: (.180446183, -.000015534, -.002290609), 1024: (.180382615, -.000005301, -.002276577),
                  2048: (.180354386, -.000001481, -.002270586), 4096: (.180347832, -.000000388, -.002269217)},
    (0.8, 0.0): {32: (0, .011451356, 0), 64: (0, .010264792, 0), 128: (0, .008500758, 0), 256: (0, .008312402, 0),
                 512: (0, .008328328, 0), 1024: (0, .008338967, 0), 2048: (0, .008341760, 0), 4096: (0, .008342211, 0)},
    (0.8, 1.25): {32: (.387842746, .027861593, -.042947499), 64: (.372921380, .025487792, -.043584427),
                  128: (.373469550, .023786371, -.043873739), 256: (.368980205, .023357651, -.042552941),
                  512: (.363747900, .023084749, -.041002228), 1024: (.360812844, .022934404, -.040136414),
                  2048: (.358281928, .022799839, -.039388829), 4096: (.357142338, .022737860, -.039051466)},
}
CFL3D = {
    (0.5, 0.0): {256: (0, .000573147, 0), 512: (0, .000131328, 0), 1024: (0, .000031615, 0), 2048: (0, .000007502, 0)},
    (0.5, 1.25): {256: (.178366720, .000649962, -.002078508), 512: (.179125033, .000163415, -.002152971),
                  1024: (.179487608, .000048451, -.002206480), 2048: (.179650543, .000020903, -.002235646)},
    (0.8, 0.0): {256: (0, .008894250, 0), 512: (0, .008479528, 0), 1024: (0, .008376058, 0), 2048: (0, .008350161, 0)},
    (0.8, 1.25): {256: (.359073197, .023411409, -.040316445), 512: (.357580694, .022902709, -.039700797),
                  1024: (.355943711, .022706750, -.039170890), 2048: (.354593186, .022606485, -.038760613)},
}
REFERENCES = {"FLO82": FLO82, "CFL3D": CFL3D}
COUNT = 1e-4  # One drag count.


def reference_case(mach, alpha_deg):
    """Key of the paper case matching the free stream, or None."""
    for key in FLO82:
        if abs(key[0] - mach) < 1e-6 and abs(key[1] - alpha_deg) < 1e-6:
            return key
    return None


def latest_study(roots=None):
    roots = roots or (Path.cwd(), Path(__file__).resolve().parent.parent)
    patterns = ("results/study_*", "build/*/results/study_*", "out/build/*/results/study_*")
    found = {path.resolve() for root in roots for pattern in patterns for path in Path(root).glob(pattern) if (path / "study.csv").is_file()}
    return max(found, key=lambda path: (path / "study.csv").stat().st_mtime, default=None)


def read_history(path):
    if not path.is_file():
        return None
    history = np.genfromtxt(path, delimiter=",", names=True)
    return np.atleast_1d(history) if history.dtype.names else None


def load_study(study, chord_ref, x_ref, exclude_radius, pressure_method):
    """One entry per mesh computed without error, ordered by increasing NC."""
    entries, failed = [], []
    with open(study / "study.csv", newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            directory = study / row["directory"]
            if row["status"] != "ok" or not (directory / "solution.vts").is_file():
                failed.append((row["directory"], row["message"] or "no solution"))
                continue
            mesh = pv.read(directory / "solution.vts")
            ni, nj, _ = mesh.dimensions
            field = lambda name, default: float(mesh.field_data[name][0]) if name in mesh.field_data else default
            history = read_history(directory / "residuals.csv")
            last = history[-1] if history is not None else None
            mach = np.asarray(mesh.cell_data["mach"], dtype=float)
            coefficients = aerodynamic_coefficients(mesh, chord_ref=chord_ref, x_ref=x_ref, pressure_method=pressure_method) or {"cl": np.nan, "cd": np.nan, "cm": np.nan}
            entries.append({
                "name": row["directory"], "mesh": mesh, "history": history,
                "nc_i": ni - 1, "nc_j": nj - 1, "cells": (ni - 1) * (nj - 1),
                "mach_inf": field("mach_inf", np.nan), "alpha_deg": np.degrees(field("alpha_rad", 0.0)),
                "iterations": int(field("iteration", last["iteration"] if last is not None else -1)),
                "rms": float(last["rms"]) if last is not None else np.nan,
                "rms_rho": float(last["rms_rho"]) if last is not None and "rms_rho" in history.dtype.names else np.nan,
                "cl": coefficients["cl"], "cd": coefficients["cd"], "cm": coefficients["cm"],
                "mach_min": float(np.nanmin(mach)), "mach_max": float(np.nanmax(mach)),
                "entropy": entropy_error(mesh),
                "entropy_excluded": entropy_error(mesh, exclude_radius) if exclude_radius > 0 else np.nan,
            })
    entries.sort(key=lambda entry: entry["cells"])
    return entries, failed


def observed_orders(entries, key):
    """Order between consecutive meshes, log(e_coarse/e_fine)/log(h_coarse/h_fine), with h = 1/sqrt(cells),
    and the least-squares slope over all meshes."""
    h = np.array([1 / np.sqrt(entry["cells"]) for entry in entries])
    e = np.array([entry[key] for entry in entries])
    local = [np.nan] + [np.log(e[n - 1] / e[n]) / np.log(h[n - 1] / h[n]) if e[n] > 0 and e[n - 1] > 0 else np.nan for n in range(1, len(e))]
    valid = (e > 0) & np.isfinite(e)
    slope = np.polyfit(np.log(h[valid]), np.log(e[valid]), 1)[0] if valid.sum() >= 2 else np.nan
    return local, slope


def write_table(entries, case, path):
    columns = ["mesh", "nc", "cells", "iterations", "rms", "rms_rho", "cl", "cd", "cm", "mach_min", "mach_max",
               "entropy_error", "entropy_order", "entropy_error_excluded", "entropy_order_excluded"]
    for name in REFERENCES:
        columns += [f"cl_{name}", f"cd_{name}", f"cm_{name}", f"dcl_{name}", f"dcd_counts_{name}", f"dcm_{name}"]
    orders, slope = observed_orders(entries, "entropy")
    orders_excluded, slope_excluded = observed_orders(entries, "entropy_excluded")
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(columns)
        for entry, order, order_excluded in zip(entries, orders, orders_excluded):
            row = [entry["name"], f"{entry['nc_i']}x{entry['nc_j']}", entry["cells"], entry["iterations"], entry["rms"], entry["rms_rho"],
                   entry["cl"], entry["cd"], entry["cm"], entry["mach_min"], entry["mach_max"],
                   entry["entropy"], order, entry["entropy_excluded"], order_excluded]
            for table in REFERENCES.values():
                values = table.get(case, {}).get(entry["nc_i"]) if case and entry["nc_i"] == entry["nc_j"] else None
                if values is None:
                    row += [""] * 6
                else:
                    row += [*values, entry["cl"] - values[0], (entry["cd"] - values[1]) / COUNT, entry["cm"] - values[2]]
            writer.writerow(row)
    return slope, slope_excluded


def print_table(entries, case, slope, slope_excluded, exclude_radius):
    print(f"{'mesh':>10} {'NC':>11} {'iter':>8} {'RES rho':>10} {'CL':>12} {'CD':>12} {'CM':>12} {'Mmax':>7} {'entropy':>10}")
    for entry in entries:
        print(f"{entry['name']:>10} {entry['nc_i']:>5}x{entry['nc_j']:<5} {entry['iterations']:>8} {entry['rms_rho']:>10.3e} "
              f"{entry['cl']:>12.9f} {entry['cd']:>12.9f} {entry['cm']:>12.9f} {entry['mach_max']:>7.4f} {entry['entropy']:>10.3e}")
        for name, table in REFERENCES.items():
            values = table.get(case, {}).get(entry["nc_i"]) if case and entry["nc_i"] == entry["nc_j"] else None
            if values is not None:
                print(f"{'':>10} {name + ' (paper)':>11} {'':>8} {'':>10} {values[0]:>12.9f} {values[1]:>12.9f} {values[2]:>12.9f}"
                      f"   dCD = {(entry['cd'] - values[1]) / COUNT:+.2f} counts")
    print(f"Entropy error: observed order (least squares over the meshes) = {slope:.3f}", end="")
    if np.isfinite(slope_excluded):
        print(f", {slope_excluded:.3f} without the cells within {exclude_radius} of the trailing edge", end="")
    print()
    if any(entry["mach_max"] >= 1 for entry in entries):
        print("Warning: supersonic pockets (Mach max >= 1). Shocks create entropy: the entropy error no longer measures the order.")
    if case is None:
        print("No paper case at this Mach number and incidence (M = 0.5 or 0.8, alpha = 0 or 1.25 deg): no reference values.")


def save_figure(figure, path, attempts=10):
    """Write to a temporary file, then rename it: on Windows an image viewer, the indexer or an
    antivirus may briefly lock a file just written (OSError), so the rename is retried."""
    temporary = path.with_name(f".{path.stem}.tmp{path.suffix}")
    figure.savefig(temporary, dpi=150)
    plt.close(figure)
    for attempt in range(attempts):
        try:
            os.replace(temporary, path)
            return
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(.2 * (attempt + 1))


def plot_residuals(entries, path):
    """Fig. 12 style (each mesh against its own iterations) and Fig. 10 style (accumulated iterations)."""
    figure, (own, accumulated) = plt.subplots(1, 2, figsize=(13, 5))
    offset = 0
    for entry in entries:
        history = entry["history"]
        if history is None:
            continue
        name = "rms_rho" if "rms_rho" in history.dtype.names else "rms"
        values = np.where(history[name] > 0, history[name], np.nan)
        label = f"{entry['nc_i']}x{entry['nc_j']}"
        own.plot(history["iteration"], np.log10(values), linewidth=1.2, label=label)
        accumulated.plot(history["iteration"] + offset, np.log10(values), linewidth=1.2, label=label)
        offset += history["iteration"][-1]
    own.set(xlabel="Iteration", ylabel="log10(RMS dRHO/dT)", title="Residual history per mesh")
    accumulated.set(xlabel="Accumulated iterations", ylabel="log10(RMS dRHO/dT)", title="Accumulated over the mesh sequence")
    for axes in (own, accumulated):
        axes.grid(True, alpha=.3)
        axes.legend()
    figure.tight_layout()
    save_figure(figure, path)


def plot_coefficients(entries, case, path):
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    nc = [entry["nc_i"] for entry in entries]
    for axis, (key, index, label) in zip(axes, (("cl", 0, "CL"), ("cd", 1, "CD"), ("cm", 2, "CM"))):
        axis.plot(nc, [entry[key] for entry in entries], "o-", label="CMFSN (JST)")
        for name, table in REFERENCES.items():
            data = table.get(case, {}) if case else {}
            if data:
                sizes = sorted(data)
                axis.plot(sizes, [data[size][index] for size in sizes], "s--", markersize=4, label=name)
        axis.set_xscale("log", base=2)
        axis.set(xlabel="NC", ylabel=label, title=label)
        axis.grid(True, alpha=.3)
        axis.legend()
    figure.tight_layout()
    save_figure(figure, path)


def plot_entropy(entries, slope, slope_excluded, exclude_radius, path):
    figure, axis = plt.subplots(figsize=(6.5, 5))
    h = np.array([1 / np.sqrt(entry["cells"]) for entry in entries])
    for key, label, order in (("entropy", "all cells", slope), ("entropy_excluded", f"r > {exclude_radius} from the TE", slope_excluded)):
        e = np.array([entry[key] for entry in entries])
        if np.any(np.isfinite(e) & (e > 0)):
            axis.loglog(h, e, "o-", label=f"{label}: order {order:.2f}")
    if len(h) >= 2 and entries[-1]["entropy"] > 0:
        guide = entries[-1]["entropy"] * (h / h[-1]) ** 2
        axis.loglog(h, guide, "k:", label="order 2")
    axis.set(xlabel="h = 1/sqrt(cells)", ylabel="L2 entropy error (s - s_inf)/s_inf", title="Spatial order from the entropy error")
    axis.grid(True, which="both", alpha=.3)
    axis.legend()
    figure.tight_layout()
    save_figure(figure, path)


def plot_solution(entry, chord_ref, pressure_method, path, mach_step=0.025):
    """Figs. 13-16 style: wall Cp (upper and lower surfaces) and Mach lines every mach_step, sonic line bold."""
    mesh = entry["mesh"]
    ni, nj, _ = mesh.dimensions
    nodes = np.asarray(mesh.points, dtype=float).reshape(nj, ni, 3)[:, :, :2]
    wall = nodes[0]
    leading, trailing = wall[np.argmin(wall[:, 0])], wall[np.argmax(wall[:, 0])]
    figure, (cp_axis, mach_axis) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.3]})

    name, x, values = wall_distribution(mesh, pressure_method)
    centres = .25 * (wall[:-1] + wall[1:] + nodes[1, :-1] + nodes[1, 1:])
    chord_line = trailing - leading
    upper = chord_line[0] * (centres[:, 1] - leading[1]) - chord_line[1] * (centres[:, 0] - leading[0]) > 0
    x_over_c = (x - leading[0]) / chord_ref
    # Hollow markers for the lower surface: both stay visible where they coincide (alpha = 0).
    cp_axis.plot(x_over_c[upper], values[upper], ".", markersize=4, label="upper surface")
    cp_axis.plot(x_over_c[~upper], values[~upper], "o", markersize=4, markerfacecolor="none", markeredgewidth=.7, label="lower surface")
    cp_axis.invert_yaxis()
    method_label = WALL_PRESSURE_LABELS[resolve_wall_pressure_method(mesh, pressure_method)]
    cp_axis.set(xlabel="x / c", ylabel="Cp" if name == "cp" else "p", title=f"Wall {'Cp' if name == 'cp' else 'pressure'} ({method_label})")
    cp_axis.grid(True, alpha=.3)
    cp_axis.legend(loc="lower right")
    text = (f"MACH {entry['mach_inf']:.3f}   ALPHA {entry['alpha_deg']:.7f}\n"
            f"CL {entry['cl']:.9f}   CD {entry['cd']:.9f}   CM {entry['cm']:.9f}\n"
            f"GRID {entry['nc_i']}x{entry['nc_j']}   NCYC {entry['iterations']}   RES {entry['rms_rho']:.3E}")
    cp_axis.text(.02, .02, text, transform=cp_axis.transAxes, family="monospace", fontsize=8, va="bottom",
                 bbox={"facecolor": "white", "alpha": .85, "edgecolor": "0.6"})

    # Mach at the cell centres; the O-mesh is closed by repeating the first column.
    mach = np.asarray(mesh.cell_data["mach"], dtype=float).reshape(nj - 1, ni - 1)
    cell_centres = .25 * (nodes[:-1, :-1] + nodes[:-1, 1:] + nodes[1:, 1:] + nodes[1:, :-1])
    close = lambda array: np.concatenate([array, array[:, :1]], axis=1)
    cx, cy, cm = close(cell_centres[..., 0]), close(cell_centres[..., 1]), close(mach)
    levels = np.round(mach_step * np.arange(1, int(entry["mach_max"] / mach_step) + 2), 10)  # Rounded: 1.0 is found exactly.
    mach_axis.contour(cx, cy, cm, levels=levels[~np.isclose(levels, 1.0)], colors="k", linewidths=.6)
    if entry["mach_min"] < 1 < entry["mach_max"]:
        mach_axis.contour(cx, cy, cm, levels=[1.0], colors="k", linewidths=2.0)
    mach_axis.fill(wall[:, 0], wall[:, 1], color="0.85", edgecolor="k", linewidth=.8)
    span = chord_ref
    mach_axis.set_xlim(leading[0] - .6 * span, leading[0] + 1.6 * span)
    mach_axis.set_ylim(-.9 * span, .9 * span)
    mach_axis.set_aspect("equal")
    mach_axis.set_title(f"Mach lines every {mach_step:g} (sonic line bold)   MIN {entry['mach_min']:.7f}   MAX {entry['mach_max']:.7f}", fontsize=9)
    figure.tight_layout()
    save_figure(figure, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("study", nargs="?", type=Path, help="results/study_NNNN (default: the most recent one)")
    parser.add_argument("--chord-ref", type=float, default=1.0, help="reference chord (paper: 1)")
    parser.add_argument("--x-ref", type=float, default=0.25, help="moment centre x (paper: 0.25)")
    parser.add_argument("--exclude-radius", type=float, default=0.05, help="entropy error also computed without the cells this close to the trailing edge (0: off)")
    parser.add_argument("--wall-pressure", default="auto", choices=WALL_PRESSURE_METHODS, help="wall pressure for Cp and the coefficients: auto = the solver's wall face pressure (rizzi for an Euler2D file with the wall_pressure_rizzi marker, cell for RANS2D or an older Euler2D file); rizzi = P_f = p1 - P_eta/2; cell = first cell; extrapolated = linear from the first two cells")
    parser.add_argument("--mach-step", type=float, default=0.025, help="spacing of the Mach lines on the solution figures (paper: 0.05)")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"), help="figure format")
    parser.add_argument("--output", type=Path, help="output directory (default: <study>/post)")
    args = parser.parse_args(argv)
    if not (args.mach_step > 0 and np.isfinite(args.mach_step)):
        parser.error("--mach-step must be positive and finite")

    study = args.study or latest_study()
    if study is None or not (study / "study.csv").is_file():
        parser.error("no study directory with a study.csv found")
    entries, failed = load_study(study, args.chord_ref, args.x_ref, args.exclude_radius, args.wall_pressure)
    for name, message in failed:
        print(f"Skipped {name}: {message}")
    if not entries:
        parser.error("no successful mesh in this study")
    cases = {(round(entry["mach_inf"], 6), round(entry["alpha_deg"], 6)) for entry in entries}
    if len(cases) != 1:
        parser.error(f"the meshes of a study must share the free stream, found {sorted(cases)}")
    case = reference_case(entries[0]["mach_inf"], entries[0]["alpha_deg"])
    output = args.output or study / "post"
    output.mkdir(parents=True, exist_ok=True)

    print(f"Study {study}  (M = {entries[0]['mach_inf']:g}, alpha = {entries[0]['alpha_deg']:g} deg, chord {args.chord_ref:g}, x_ref {args.x_ref:g}, wall pressure {args.wall_pressure})")
    slope, slope_excluded = write_table(entries, case, output / "table.csv")
    print_table(entries, case, slope, slope_excluded, args.exclude_radius)
    plot_residuals(entries, output / f"residuals.{args.format}")
    plot_coefficients(entries, case, output / f"coefficients.{args.format}")
    plot_entropy(entries, slope, slope_excluded, args.exclude_radius, output / f"entropy.{args.format}")
    for entry in entries:
        plot_solution(entry, args.chord_ref, args.wall_pressure, output / f"solution_{entry['name']}.{args.format}", args.mach_step)
    print(f"Written to {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
