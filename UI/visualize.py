"""Inspect CMFSN structured-grid exports without interpolating cell values.

Install: python -m pip install -r UI/requirements.txt
Latest:  python UI/visualize.py                 (most recent results/run_*/solution)
Series:  python UI/visualize.py results/run_0001/solution.pvd
Single:  python UI/visualize.py results/run_0001/solution.vts
Zoom:    python UI/visualize.py results/run_0001/solution.vts --bounds -0.25 1.25 -0.75 0.75
Image:   python UI/visualize.py results/run_0001/solution.vts --screenshot pressure.png
Keys:    Space play/pause, F5 reload, Ctrl+O open, Ctrl+S export image
Theme:   python UI/visualize.py --theme clair    (or the button in the header; the choice is remembered)
Export:  python UI/visualize.py results/run_0001/solution.pvd --gif pression.gif --fps 15
         python UI/visualize.py results/run_0001/solution.pvd --curve coefficients --save-curve coefficients.png   (.svg, .pdf, .csv)
Hover:   the flow view gives the values of the cell under the cursor; a curve gives its nearest point
"""
from __future__ import annotations

import argparse
from functools import partial
from math import ceil
from time import perf_counter
from pathlib import Path
from string import Template
import sys

import numpy as np
import pyvista as pv


FIELDS = {
    "pressure": "Pression",
    "mach": "Mach",
    "density": "Densité",
    "temperature": "Température (p/ρ)",
    "cp": "Coefficient de pression Cp",
    "velocity_u": "Vitesse u",
    "velocity_v": "Vitesse v",
    "speed": "Norme de la vitesse",
    "energy": "Énergie totale / volume",
}


# Every colour of the interface, the flow view and the curves comes from the active theme.
THEMES = {
    "sombre": {
        "name": "sombre", "window": "#0a111e", "card": "#121b2d", "card_border": "#1f2c44", "tile": "#0d1626",
        "header_start": "#0d1a30", "header_end": "#0f4a52", "text": "#dfe7f3", "title": "#f4f7fc", "muted": "#8b9cb6", "section": "#62779a",
        "accent": "#2dd4bf", "accent_hover": "#5eead4", "accent_text": "#04211d", "input": "#0d1626", "input_border": "#283854",
        "button": "#17233a", "button_hover": "#1f2e4b", "track": "#283854", "badge_bg": "#0f332f", "badge_text": "#5eead4",
        "plot": "#0d1524", "grid": "#1d2a42", "spine": "#2a3a57", "curve": "#2dd4bf", "marker": "#fbbf24", "nan": "#334155", "edge": "#5b6b85",
    },
    "clair": {
        "name": "clair", "window": "#eef2f7", "card": "#ffffff", "card_border": "#e1e8f1", "tile": "#f5f8fb",
        "header_start": "#0f1f3a", "header_end": "#136068", "text": "#1f2a3d", "title": "#14243c", "muted": "#6b7a90", "section": "#8090a6",
        "accent": "#0f9488", "accent_hover": "#0c7d73", "accent_text": "#ffffff", "input": "#ffffff", "input_border": "#d6dfeb",
        "button": "#f5f8fb", "button_hover": "#e9eff6", "track": "#dbe4ee", "badge_bg": "#e6f6f3", "badge_text": "#0f766e",
        "plot": "#ffffff", "grid": "#e8edf4", "spine": "#dce4ef", "curve": "#0f9488", "marker": "#d98b32", "nan": "#cbd5e1", "edge": "#66748a",
    },
}
THEME = dict(THEMES["sombre"])  # Active theme, updated in place by configure_theme.

STYLESHEET = Template("""
    QWidget { font-family: 'Segoe UI', 'DejaVu Sans'; font-size: 13px; color: $text; }
    QMainWindow, #root { background: $window; }
    #header { background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 $header_start, stop:1 $header_end); }
    #brand { color: white; font-size: 25px; font-weight: 700; }
    #subtitle { color: #9fb3cc; font-size: 11px; }
    #filename { color: #e2ebf5; background: rgba(255, 255, 255, 0.10); border-radius: 12px; padding: 5px 14px; }
    #sidebar, #flowCard, #curveCard { background: $card; border: 1px solid $card_border; border-radius: 14px; }
    #section { color: $section; font-size: 10px; font-weight: 700; }
    #title { font-size: 19px; font-weight: 600; color: $title; }
    #tile { background: $tile; border: 1px solid $card_border; border-radius: 9px; }
    #metric { font-size: 16px; font-weight: 700; color: $title; }
    #executionMetric { font-size: 18px; font-weight: 600; color: $title; }
    #coefficients { font-size: 14px; font-weight: 600; color: $title; background: $tile; border: 1px solid $card_border; border-radius: 9px; padding: 5px 12px; }
    #muted { color: $muted; font-size: 12px; }
    #badge { background: $badge_bg; color: $badge_text; border-radius: 6px; padding: 6px 10px; font-size: 10px; font-weight: 600; }
    QPushButton { border: 1px solid $input_border; border-radius: 8px; background: $button; padding: 9px 12px; font-weight: 600; }
    QPushButton:hover { background: $button_hover; border-color: $accent; }
    QPushButton#primary { background: $accent; color: $accent_text; border: none; }
    QPushButton#primary:hover { background: $accent_hover; }
    QPushButton#ghost { background: rgba(255, 255, 255, 0.08); color: white; border: 1px solid rgba(255, 255, 255, 0.22); }
    QPushButton#ghost:hover { background: rgba(255, 255, 255, 0.18); border-color: rgba(255, 255, 255, 0.45); }
    QSpinBox, QDoubleSpinBox { border: 1px solid $input_border; border-radius: 6px; background: $input; padding: 5px 8px; }
    QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { width: 0px; border: none; }
    QComboBox { border: 1px solid $input_border; border-radius: 8px; background: $input; padding: 9px; }
    QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: $accent; }
    QComboBox::drop-down { width: 24px; border: none; }
    QCheckBox { spacing: 8px; padding: 6px 0; }
    QCheckBox::indicator { width: 15px; height: 15px; border-radius: 4px; border: 1px solid $input_border; background: $input; }
    QCheckBox::indicator:hover { border-color: $accent; }
    QCheckBox::indicator:checked { background: $accent; border-color: $accent; }
    QWidget:disabled { color: $muted; }
    QPushButton:disabled { background: $tile; border-color: $card_border; color: $muted; }
    QComboBox QAbstractItemView { background: $card; color: $text; border: 1px solid $card_border; selection-background-color: $accent; selection-color: $accent_text; }
    QToolTip { background: $card; color: $text; border: 1px solid $input_border; padding: 6px; }
    QSlider::groove:horizontal { background: $track; height: 6px; border-radius: 3px; }
    QSlider::sub-page:horizontal { background: $accent; border-radius: 3px; }
    QSlider::handle:horizontal { background: $accent; border: 2px solid $card; width: 14px; margin: -6px 0; border-radius: 8px; }
    QSlider::handle:horizontal:disabled { background: $muted; }
""")

# Short names used when reading the values of the cell under the cursor.
PROBE_LABELS = {"pressure": "p", "mach": "M", "density": "ρ", "temperature": "θ", "cp": "Cp", "velocity_u": "u", "velocity_v": "v"}
COEFFICIENTS = (("cl", "Cl"), ("cd", "Cd"), ("cm", "Cm"))
MAX_COEFFICIENT_SAMPLES = 200  # States read from a series to draw the evolution of the coefficients.
MAX_GIF_FRAMES = 300  # A longer series is sampled evenly.
CURVES = ("convergence", "paroi", "coefficients")  # Order of the curve selector.


def gif_frames(count, maximum=MAX_GIF_FRAMES, stride=1):
    """Indices of the states written to a GIF: every `stride`-th state, then evenly thinned to `maximum`; the last state is kept."""
    indices = np.arange(0, count, max(1, int(stride)))
    if indices.size and indices[-1] != count - 1:
        indices = np.append(indices, count - 1)
    if indices.size > maximum:
        indices = indices[np.unique(np.linspace(0, indices.size - 1, maximum).round().astype(int))]
    return indices


def write_gif(path, images, fps):
    """Write RGB(A) arrays as a looping GIF (Pillow ships with matplotlib)."""
    from PIL import Image
    if not images:
        raise ValueError("Aucune image à écrire")
    frames = [Image.fromarray(np.asarray(image)[..., :3]).quantize(colors=256, method=Image.Quantize.MEDIANCUT) for image in images]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=max(20, round(1000 / fps)), loop=0, optimize=False, disposal=2)


def prepare_solution(mesh) -> pv.StructuredGrid:
    if not isinstance(mesh, pv.StructuredGrid) or mesh.dimensions[2] != 1:
        raise ValueError("Le fichier doit contenir une grille structuree 2D (.vts).")
    required = {"pressure", "density", "mach", "energy", "velocity"}
    if not required.issubset(mesh.cell_data):
        raise ValueError(f"Champs manquants : {', '.join(sorted(required - set(mesh.cell_data)))}")
    velocity = np.asarray(mesh.cell_data["velocity"])
    if velocity.shape != (mesh.n_cells, 3):
        raise ValueError("Le champ velocity doit avoir trois composantes par cellule.")
    mesh.cell_data["velocity_u"] = velocity[:, 0]
    mesh.cell_data["velocity_v"] = velocity[:, 1]
    mesh.cell_data["speed"] = np.linalg.norm(velocity, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        mesh.cell_data["temperature"] = np.asarray(mesh.cell_data["pressure"]) / np.asarray(mesh.cell_data["density"])
    return mesh


def latest_solution(roots=None):
    """Most recently written run, searched from the working directory and the repository."""
    roots = roots or (Path.cwd(), Path(__file__).resolve().parent.parent)
    patterns = ("results/run_*/solution.*", "build/*/results/run_*/solution.*", "out/build/*/results/run_*/solution.*")
    found = {path.resolve() for root in roots for pattern in patterns for path in Path(root).glob(pattern) if path.suffix.lower() in (".pvd", ".vts")}
    return max(found, key=lambda path: path.stat().st_mtime, default=None)


def flow_summary(mesh) -> str:
    """Free-stream conditions stored with the solution, when present."""
    parts = []
    if "mach_inf" in mesh.field_data:
        parts.append(f"M∞ = {float(mesh.field_data['mach_inf'][0]):.4g}")
    if "alpha_rad" in mesh.field_data:
        parts.append(f"α = {np.degrees(float(mesh.field_data['alpha_rad'][0])):.4g}°")
    return " · ".join(parts)


def aerodynamic_coefficients(mesh):
    """Force coefficients on the wall (j = 0), integrated from the row of cells touching it.

    Pressure and velocity are taken at the first cell centres. The friction part is present only
    when the solution stores a Reynolds number: tau = mu * tangential velocity / wall distance,
    with the solver's mu = M_inf * sqrt(gamma) / Re. Reference length: the chord (x extent of the
    wall). Cm is about the quarter chord, positive nose-up. Returns None for a free stream at rest.
    """
    data = mesh.field_data
    value = lambda name, default: float(data[name][0]) if name in data else default
    mach, gamma, alpha = value("mach_inf", 0.0), value("gamma", 1.4), value("alpha_rad", 0.0)
    rho_inf, p_inf, reynolds = value("rho_inf", 1.0), value("pressure_inf", 1.0), value("reynolds_inf", np.inf)
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
    pressure_force = (row("pressure")[:, 0] - p_inf)[:, None] * normal
    friction_force = np.zeros_like(pressure_force)
    viscous = bool(np.isfinite(reynolds) and reynolds > 0)
    if viscous:
        distance = np.abs(np.einsum("ij,ij->i", normal, middle - centers))  # Wall distance times face length.
        valid = distance > 0
        unit = normal[valid] / length[valid, None]
        velocity = row("velocity")[valid, :2]
        tangential = velocity - np.einsum("ij,ij->i", velocity, unit)[:, None] * unit
        viscosity = rho_inf * mach * np.sqrt(gamma * p_inf / rho_inf) / reynolds
        friction_force[valid] = viscosity * tangential * (length[valid] ** 2 / distance[valid])[:, None]
    wall_x = points[0, :, 0]
    chord = float(np.ptp(wall_x))
    if not chord > 0:
        return None
    origin = np.array([wall_x.min() + .25 * chord, points[0, np.argmin(wall_x), 1]])
    lever = middle - origin

    def coefficients(force):
        fx, fy = np.nansum(force, axis=0) / (dynamic_pressure * chord)
        moment = -np.nansum(lever[:, 0] * force[:, 1] - lever[:, 1] * force[:, 0]) / (dynamic_pressure * chord * chord)
        return np.array([-fx * np.sin(alpha) + fy * np.cos(alpha), fx * np.cos(alpha) + fy * np.sin(alpha), moment])

    pressure_part, friction_part = coefficients(pressure_force), coefficients(friction_force)
    total = pressure_part + friction_part
    return {"cl": total[0], "cd": total[1], "cm": total[2], "cd_pressure": pressure_part[1], "cd_friction": friction_part[1], "viscous": viscous, "chord": chord}


def wall_distribution(mesh):
    """Cp (pressure when Cp is absent) in the row of cells touching the wall (j = 0), against x."""
    name = "cp" if "cp" in mesh.cell_data else "pressure"
    ni, nj, _ = mesh.dimensions
    x = np.asarray(mesh.points).reshape(nj, ni, 3)[:2, :, 0]
    centers = .25 * (x[0, :-1] + x[0, 1:] + x[1, :-1] + x[1, 1:])
    return name, centers, np.asarray(mesh.cell_data[name]).reshape(nj - 1, ni - 1)[0]


def load_solution(path: Path) -> pv.StructuredGrid:
    return prepare_solution(pv.read(path))


class SolutionSequence:
    """Keep only the time index; load the selected physical state on demand."""
    def __init__(self, path: Path):
        self.path = Path(path).resolve()
        self.reader = None
        if self.path.suffix.lower() == ".pvd":
            self.reader = pv.PVDReader(self.path)
            self.times = np.asarray(self.reader.time_values, dtype=float)
            if self.times.size == 0 or not np.all(np.isfinite(self.times)) or np.any(np.diff(self.times) <= 0):
                raise ValueError("La série doit contenir des temps physiques finis et strictement croissants.")
        elif self.path.suffix.lower() == ".vts":
            self.times = np.array([np.nan])
        else:
            raise ValueError("Choisir une solution .vts ou une série temporelle .pvd.")

    def __len__(self):
        return self.times.size

    def read(self, index):
        if not 0 <= index < len(self):
            raise ValueError("Indice temporel hors de la série.")
        if self.reader is None:
            mesh = load_solution(self.path)
            if "TimeValue" in mesh.field_data:
                self.times[0] = float(mesh.field_data["TimeValue"][0])
            return mesh
        self.reader.set_active_time_point(index)
        blocks = self.reader.read()
        if len(blocks) != 1:
            raise ValueError("Chaque instant doit contenir une seule grille structurée.")
        return prepare_solution(blocks[0])


def load_history(path: Path) -> np.ndarray:
    history = np.genfromtxt(path, delimiter=",", names=True, ndmin=1, encoding="utf-8")
    if history.dtype.names is None or not {"iteration", "rms"}.issubset(history.dtype.names):
        raise ValueError("L'historique doit contenir les colonnes iteration et rms.")
    if history.size == 0 or not np.all(np.isfinite(history["iteration"])) or not np.all(np.isfinite(history["rms"])) or np.any(history["rms"] < 0):
        raise ValueError("Historique vide ou contenant un residu invalide.")
    return history


def execution_time_text(history) -> str:
    """Format the last recorded wall-clock duration; never infer it from physical time."""
    if history is None or "elapsed_seconds" not in (history.dtype.names or ()) or history.size == 0:
        return "Indisponible"
    elapsed = float(history["elapsed_seconds"][-1])
    if not np.isfinite(elapsed) or elapsed < 0:
        return "Indisponible"
    total_ms = round(elapsed * 1000)
    hours, remainder = divmod(total_ms, 3_600_000)
    minutes, remainder = divmod(remainder, 60_000)
    seconds, milliseconds = divmod(remainder, 1000)
    clock = f"{seconds:02d},{milliseconds:03d} s"
    if hours:
        return f"{hours} h {minutes:02d} min {clock}"
    if minutes:
        return f"{minutes} min {clock}"
    return f"{seconds},{milliseconds:03d} s"


class SolutionViewer:
    def __init__(self, mesh, *, field="pressure", bounds=None, off_screen=False, plotter=None):
        if field not in FIELDS or field not in mesh.cell_data:
            raise ValueError(f"Champ indisponible : {field}")
        self.mesh = mesh
        self.vectors_enabled = False
        self.vector_maximum = 600
        self.vector_scale = 1.0
        self.vector_count = 0
        self.vector_actor = None
        self.vector_cell_ids = np.empty(0, dtype=np.int64)
        self._vector_centers = None
        self._vector_sample_key = None
        self.fields = [name for name in FIELDS if name in mesh.cell_data]
        self.plotter = plotter if plotter is not None else pv.Plotter(off_screen=off_screen, window_size=(1100, 570))
        self.plotter.clear()
        self.actor = self.plotter.add_mesh(mesh, scalars=field, preference="cell", cmap="viridis", lighting=False, show_edges=False, nan_color=THEME["nan"], scalar_bar_args={"title": "", "vertical": False, "position_x": .22, "position_y": .025, "width": .56, "height": .10, "color": THEME["text"], "label_font_size": 12, "n_labels": 5})
        self.apply_theme()
        self.plotter.view_xy()
        self.plotter.enable_parallel_projection()
        self.plotter.enable_2d_style()
        if bounds is not None:
            self.focus(bounds)
        self.set_field(field)

    def set_mesh(self, mesh, *, field=None):
        """Replace a time sample without rebuilding the scene or its camera."""
        field = field or self.current_field
        if field not in FIELDS or field not in mesh.cell_data:
            raise ValueError(f"Champ indisponible : {field}")
        if not np.any(np.isfinite(mesh.cell_data[field])):
            raise ValueError(f"Aucune valeur finie pour {field}")
        self.mesh = mesh
        self.fields = [name for name in FIELDS if name in mesh.cell_data]
        self.actor.mapper.dataset = mesh
        self._vector_centers = None
        self._vector_sample_key = None
        self.set_field(field)

    def apply_theme(self):
        """Colours of the view that follow the interface theme: background, colour bar, edges, missing values."""
        rgb = lambda key: pv.Color(THEME[key]).float_rgb
        self.plotter.set_background(THEME["plot"])
        scalar_bar = self.plotter.scalar_bar
        scalar_bar.GetBackgroundProperty().SetColor(*rgb("plot"))
        scalar_bar.SetDrawBackground(True)
        labels = scalar_bar.GetLabelTextProperty()
        labels.SetColor(*rgb("text"))
        labels.SetBackgroundColor(*rgb("plot"))
        labels.SetBackgroundOpacity(1.0)
        self.actor.mapper.lookup_table.nan_color = THEME["nan"]
        self.actor.prop.edge_color = THEME["edge"]
        self.plotter.render()

    def find_cell(self, x, y):
        """Index of the cell containing (x, y), or -1 outside the mesh. The locator is built once per state."""
        if getattr(self, "_locator_mesh", None) is not self.mesh:
            from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator
            self._locator = vtkStaticCellLocator()
            self._locator.SetDataSet(self.mesh)
            self._locator.BuildLocator()
            self._locator_mesh = self.mesh
        return int(self._locator.FindCell((float(x), float(y), 0.0)))

    def probe(self, x, y):
        """Values of the cell under (x, y) as one line of text, or None outside the mesh."""
        cell = self.find_cell(x, y) if np.isfinite(x) and np.isfinite(y) else -1
        if cell < 0:
            return None
        columns = self.mesh.dimensions[0] - 1
        parts = [f"x = {x:.5g}, y = {y:.5g}", f"cellule ({cell % columns}, {cell // columns})"]
        parts += [f"{label} = {float(self.mesh.cell_data[name][cell]):.6g}" for name, label in PROBE_LABELS.items() if name in self.mesh.cell_data]
        return "  ·  ".join(parts)

    def set_field(self, field):
        values = np.asarray(self.mesh.cell_data[field])
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            raise ValueError(f"Aucune valeur finie pour {field}")
        low, high = float(finite.min()), float(finite.max())
        padding = max(abs(low), 1)*1e-6
        self.mesh.set_active_scalars(field, preference="cell")
        self.actor.mapper.array_name = field
        self.actor.mapper.scalar_range = (low, high) if low < high else (low-padding, high+padding)
        self.plotter.scalar_bar.SetLookupTable(self.actor.mapper.lookup_table)
        self.current_field = field
        self.stats = {"min": low, "max": high, "invalid": values.size-finite.size}
        self.plotter.render()

    def set_edges(self, visible):
        self.actor.prop.show_edges = bool(visible)
        self.actor.prop.edge_color = THEME["edge"]
        self.actor.prop.line_width = .5
        self.plotter.render()

    def set_colormap(self, name):
        self.actor.mapper.lookup_table.apply_cmap(name)
        self.plotter.render()

    def reset_camera(self):
        self.plotter.view_xy()
        self.plotter.reset_camera()
        self.update_velocity_vectors()

    def focus(self, bounds):
        xmin, xmax, ymin, ymax = bounds
        if xmin >= xmax or ymin >= ymax or not np.all(np.isfinite(bounds)):
            raise ValueError("Bornes invalides : xmin < xmax et ymin < ymax requis")
        self.plotter.reset_camera(bounds=(xmin, xmax, ymin, ymax, 0, 0))
        self.update_velocity_vectors()


    def _visible_xy_bounds(self):
        height = 2 * float(self.plotter.camera.parallel_scale)
        width_pixels, height_pixels = self.plotter.window_size
        width = height * max(width_pixels, 1) / max(height_pixels, 1)
        x, y, _ = self.plotter.camera.focal_point
        return (x-width/2, x+width/2, y-height/2, y+height/2)

    def _sample_velocity(self, bounds, maximum):
        if self._vector_centers is None:
            ni, nj, _ = self.mesh.dimensions
            points = np.asarray(self.mesh.points).reshape(nj, ni, 3)
            self._vector_centers = ((points[:-1, :-1, :2] + points[:-1, 1:, :2] + points[1:, :-1, :2] + points[1:, 1:, :2]) * .25).reshape(-1, 2)
            velocity = np.asarray(self.mesh.cell_data["velocity"])
            speeds = np.asarray(self.mesh.cell_data["speed"])
            self._vector_valid = np.isfinite(self._vector_centers).all(axis=1) & np.isfinite(velocity).all(axis=1) & np.isfinite(speeds) & (speeds > 0)
            self._vector_reference = float(np.max(speeds[self._vector_valid])) if np.any(self._vector_valid) else 0.0
        key = tuple(bounds) + (maximum,)
        if key == self._vector_sample_key:
            return self.vector_cell_ids, self._vector_spacing
        xmin, xmax, ymin, ymax = bounds
        width, height = xmax-xmin, ymax-ymin
        if width <= 0 or height <= 0 or not np.all(np.isfinite(bounds)):
            return np.empty(0, dtype=np.int64), 0.0
        centers = self._vector_centers
        inside = self._vector_valid & (centers[:, 0] >= xmin) & (centers[:, 0] < xmax) & (centers[:, 1] >= ymin) & (centers[:, 1] < ymax)
        ids = np.flatnonzero(inside)
        nx = min(maximum, max(1, int(np.sqrt(maximum * width / height))))
        ny = max(1, maximum // nx)
        if ids.size:
            x_bin = np.minimum(nx-1, ((centers[ids, 0]-xmin) * nx / width).astype(np.int64))
            y_bin = np.minimum(ny-1, ((centers[ids, 1]-ymin) * ny / height).astype(np.int64))
            _, first = np.unique(x_bin + nx*y_bin, return_index=True)
            ids = ids[first]
        self._vector_sample_key = key
        self.vector_cell_ids = ids
        self._vector_spacing = min(width/nx, height/ny)
        return ids, self._vector_spacing

    def set_velocity_vectors(self, enabled, *, maximum=600, scale=1.0):
        if maximum < 1 or not np.isfinite(scale) or scale <= 0:
            raise ValueError("Le nombre de flèches et leur échelle doivent être positifs")
        self.vectors_enabled = bool(enabled)
        self.vector_maximum = int(maximum)
        self.vector_scale = float(scale)
        self.update_velocity_vectors()

    def update_velocity_vectors(self):
        self.vector_count = 0
        if not self.vectors_enabled:
            self.plotter.remove_actor("velocity_vectors", reset_camera=False, render=False)
            self.vector_actor = None
            self.plotter.render()
            return
        bounds = self._visible_xy_bounds()
        ids, spacing = self._sample_velocity(bounds, self.vector_maximum)
        if ids.size == 0 or self._vector_reference <= 0:
            self.plotter.remove_actor("velocity_vectors", reset_camera=False, render=False)
            self.vector_actor = None
            self.plotter.render()
            return
        anchors = np.zeros((ids.size, 3), dtype=self._vector_centers.dtype)
        anchors[:, :2] = self._vector_centers[ids]
        anchors[:, 2] = max(bounds[1]-bounds[0], bounds[3]-bounds[2]) * 1e-4
        cloud = pv.PolyData(anchors)
        cloud.point_data["velocity"] = np.asarray(self.mesh.cell_data["velocity"])[ids]
        cloud.point_data["speed"] = np.asarray(self.mesh.cell_data["speed"])[ids]
        arrow_points = np.array([[0, -.045, 0], [.65, -.045, 0], [.65, -.14, 0], [1, 0, 0], [.65, .14, 0], [.65, .045, 0], [0, .045, 0]])
        arrow = pv.PolyData(arrow_points, faces=np.array([7, 0, 1, 2, 3, 4, 5, 6]))
        self.vector_length_factor = .7 * spacing * self.vector_scale / self._vector_reference
        self.vector_glyphs = cloud.glyph(orient="velocity", scale="speed", factor=self.vector_length_factor, geom=arrow)
        if self.vector_actor is None:
            self.vector_actor = self.plotter.add_mesh(self.vector_glyphs, name="velocity_vectors", color="#ffffff", show_edges=True, edge_color="#18344e", line_width=.7, lighting=False, show_scalar_bar=False, reset_camera=False, render=False)
        else:
            self.vector_actor.mapper.dataset = self.vector_glyphs
        self.vector_count = int(ids.size)
        self.plotter.render()


def build_stylesheet() -> str:
    return STYLESHEET.substitute(THEME)


def saved_theme() -> str:
    from PySide6 import QtCore
    name = QtCore.QSettings("CMFSN", "visualize").value("theme", "sombre")
    return name if name in THEMES else "sombre"


def configure_theme(app, name=None):
    """Activate a theme (the current one by default) with an explicit palette, instead of inheriting the desktop's."""
    from PySide6 import QtCore, QtGui
    if name is not None:
        THEME.update(THEMES[name])
    app.setStyle("Fusion")
    app.styleHints().setColorScheme(QtCore.Qt.ColorScheme.Dark if THEME["name"] == "sombre" else QtCore.Qt.ColorScheme.Light)
    palette = QtGui.QPalette()
    role = QtGui.QPalette.ColorRole
    colors = {
        role.Window: "window", role.WindowText: "text", role.Base: "input", role.AlternateBase: "tile",
        role.ToolTipBase: "card", role.ToolTipText: "text", role.Text: "text", role.Button: "button", role.ButtonText: "text",
        role.BrightText: "title", role.Highlight: "accent", role.HighlightedText: "accent_text", role.PlaceholderText: "muted", role.Accent: "accent",
    }
    for key, color in colors.items():
        palette.setColor(key, QtGui.QColor(THEME[color]))
    for key in [role.WindowText, role.Text, role.ButtonText, role.PlaceholderText]:
        palette.setColor(QtGui.QPalette.ColorGroup.Disabled, key, QtGui.QColor(THEME["muted"]))
    palette.setColor(QtGui.QPalette.ColorGroup.Disabled, role.Base, QtGui.QColor(THEME["tile"]))
    app.setPalette(palette)


def create_window(solution, *, field="pressure", history_path=None, bounds=None, headless=False, vectors=False, frame=0):
    import os
    os.environ.setdefault("QT_API", "pyside6")
    from PySide6 import QtCore, QtGui, QtWidgets
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    from matplotlib import get_data_path

    configure_theme(QtWidgets.QApplication.instance())

    # Bundle-readable fonts also work on Windows without system-font discovery.
    for filename in ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf"):
        QtGui.QFontDatabase.addApplicationFont(str(Path(get_data_path()) / "fonts" / "ttf" / filename))

    class ViewerWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self.resize(1440, 980)
            self.setWindowTitle("CMFSN · Visualisation des résultats")
            self.setStyleSheet(build_stylesheet())
            root = QtWidgets.QWidget(objectName="root")
            outer = QtWidgets.QVBoxLayout(root)
            outer.setContentsMargins(0, 0, 0, 0)
            outer.setSpacing(0)
            self.setCentralWidget(root)
            header = QtWidgets.QFrame(objectName="header")
            top = QtWidgets.QHBoxLayout(header)
            top.setContentsMargins(24, 16, 24, 16)
            top.setSpacing(10)
            brand = QtWidgets.QVBoxLayout()
            brand.setSpacing(0)
            brand.addWidget(QtWidgets.QLabel("CMFSN", objectName="brand"))
            brand.addWidget(QtWidgets.QLabel("POST-TRAITEMENT  /  ÉCOULEMENT 2D", objectName="subtitle"))
            top.addLayout(brand)
            top.addStretch()
            self.file_label = QtWidgets.QLabel(objectName="filename")
            top.addWidget(self.file_label)
            self.theme_button = QtWidgets.QPushButton(objectName="ghost")
            self.theme_button.setToolTip("Basculer entre le thème sombre et le thème clair. Le choix est conservé.")
            self.theme_button.clicked.connect(lambda: self.apply_theme("clair" if THEME["name"] == "sombre" else "sombre", save=True))
            top.addWidget(self.theme_button)
            reload_button = QtWidgets.QPushButton("Actualiser", objectName="ghost")
            reload_button.setToolTip("Relire le fichier et l'historique (F5), par exemple pendant un calcul en cours. La vue est conservée.")
            reload_button.clicked.connect(self.reload)
            top.addWidget(reload_button)
            open_button = QtWidgets.QPushButton("Ouvrir une solution…", objectName="primary")
            open_button.clicked.connect(self.open_file)
            top.addWidget(open_button)
            outer.addWidget(header)
            body = QtWidgets.QHBoxLayout()
            body.setContentsMargins(20, 20, 20, 12)
            body.setSpacing(18)
            outer.addLayout(body, 1)
            sidebar = QtWidgets.QFrame(objectName="sidebar")
            sidebar.setFixedWidth(250)
            side = QtWidgets.QVBoxLayout(sidebar)
            side.setContentsMargins(20, 22, 20, 20)
            side.setSpacing(12)
            side.addWidget(QtWidgets.QLabel("AFFICHAGE", objectName="section"))
            side.addWidget(QtWidgets.QLabel("Grandeur physique"))
            self.field_combo = QtWidgets.QComboBox()
            self.field_combo.currentIndexChanged.connect(self.change_field)
            side.addWidget(self.field_combo)
            side.addWidget(QtWidgets.QLabel("Palette de couleurs"))
            self.palette_combo = QtWidgets.QComboBox()
            self.palette_combo.addItems(["viridis", "turbo", "plasma", "magma", "coolwarm", "cividis"])
            self.palette_combo.currentTextChanged.connect(self.change_palette)
            side.addWidget(self.palette_combo)
            self.edges = QtWidgets.QCheckBox("Afficher les arêtes")
            self.edges.toggled.connect(self.change_edges)
            side.addWidget(self.edges)
            fit = QtWidgets.QPushButton("Domaine entier")
            fit.clicked.connect(self.fit_domain)
            side.addWidget(fit)
            focus = QtWidgets.QPushButton("Zoom profil NACA")
            focus.setToolTip("Vue autour d'un profil normalisé de corde 0 à 1.")
            focus.clicked.connect(self.focus_airfoil)
            side.addWidget(focus)
            self.color_lock = QtWidgets.QCheckBox("Fixer l'échelle de couleurs")
            self.color_lock.setToolTip("Conserver la même échelle entre les instants pour le champ choisi.")
            self.color_lock.toggled.connect(self.change_color_lock)
            side.addWidget(self.color_lock)
            side.addSpacing(12)
            side.addWidget(QtWidgets.QLabel("VALEURS DU CHAMP", objectName="section"))
            self.min_label = QtWidgets.QLabel(objectName="metric")
            self.max_label = QtWidgets.QLabel(objectName="metric")
            for caption, label in (("Minimum", self.min_label), ("Maximum", self.max_label)):
                tile = QtWidgets.QFrame(objectName="tile")
                box = QtWidgets.QHBoxLayout(tile)
                box.setContentsMargins(12, 9, 12, 9)
                box.addWidget(QtWidgets.QLabel(caption, objectName="muted"))
                box.addStretch()
                box.addWidget(label)
                side.addWidget(tile)
            self.validity = QtWidgets.QLabel(objectName="muted")
            self.validity.setWordWrap(True)
            side.addWidget(self.validity)
            side.addSpacing(12)
            side.addWidget(QtWidgets.QLabel("TEMPS D’EXÉCUTION", objectName="section"))
            self.execution_time_label = QtWidgets.QLabel(objectName="executionMetric")
            self.execution_time_label.setWordWrap(True)
            self.execution_time_label.setToolTip("Durée réelle cumulée du calcul, à la dernière mesure de l’historique. Elle reste fixe quand on change l’instant physique affiché.")
            side.addWidget(self.execution_time_label)
            self.execution_time_caption = QtWidgets.QLabel(objectName="muted")
            self.execution_time_caption.setWordWrap(True)
            side.addWidget(self.execution_time_caption)
            side.addStretch()
            self.mesh_label = QtWidgets.QLabel(objectName="muted")
            self.mesh_label.setWordWrap(True)
            side.addWidget(self.mesh_label)
            image_button = QtWidgets.QPushButton("Exporter une image", objectName="primary")
            image_button.setToolTip("Capture de toute la fenêtre (Ctrl+S).")
            image_button.clicked.connect(self.save_image)
            side.addWidget(image_button)
            self.gif_button = QtWidgets.QPushButton("Exporter un GIF")
            self.gif_button.setToolTip("Animation de la vue du champ sur les états de la série (calcul instationnaire), avec une échelle de couleurs commune.")
            self.gif_button.clicked.connect(self.save_gif)
            side.addWidget(self.gif_button)
            body.addWidget(sidebar)
            center = QtWidgets.QVBoxLayout()
            center.setSpacing(16)
            body.addLayout(center, 1)
            flow = QtWidgets.QFrame(objectName="flowCard")
            layout = QtWidgets.QVBoxLayout(flow)
            layout.setContentsMargins(16, 12, 16, 8)
            row = QtWidgets.QHBoxLayout()
            self.field_title = QtWidgets.QLabel(objectName="title")
            row.addWidget(self.field_title)
            row.addStretch()
            tip = "Efforts sur la paroi (j = 0) intégrés depuis la première rangée de cellules, pour l'état affiché.\nLongueur de référence : la corde. Cm au quart de corde, positif à cabrer.\nLe frottement est estimé par μ · vitesse tangentielle / distance à la paroi."
            self.coefficients_label = QtWidgets.QLabel(objectName="coefficients")
            self.coefficients_label.setToolTip(tip)
            row.addWidget(self.coefficients_label)
            self.coefficients_caption = QtWidgets.QLabel(objectName="muted")
            self.coefficients_caption.setToolTip(tip)
            row.addWidget(self.coefficients_caption)
            row.addSpacing(12)
            row.addWidget(QtWidgets.QLabel("VALEURS AUX CELLULES", objectName="badge"))
            layout.addLayout(row)
            timeline = QtWidgets.QHBoxLayout()
            self.play_button = QtWidgets.QPushButton("Lire")
            self.play_button.setFixedWidth(80)
            self.play_button.clicked.connect(self.toggle_playback)
            timeline.addWidget(self.play_button)
            previous = QtWidgets.QPushButton("‹")
            previous.setFixedWidth(36)
            previous.clicked.connect(lambda: self.frame_slider.setValue(max(0, self.frame_index - 1)))
            timeline.addWidget(previous)
            self.frame_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
            self.frame_slider.setTracking(False)
            self.frame_slider.valueChanged.connect(self.select_frame)
            timeline.addWidget(self.frame_slider, 1)
            following = QtWidgets.QPushButton("›")
            following.setFixedWidth(36)
            following.clicked.connect(lambda: self.frame_slider.setValue(min(len(self.sequence) - 1, self.frame_index + 1)))
            timeline.addWidget(following)
            self.frame_spin = QtWidgets.QSpinBox()
            self.frame_spin.setPrefix("État ")
            self.frame_spin.setFixedWidth(105)
            self.frame_spin.valueChanged.connect(lambda value: self.frame_slider.setValue(value - 1))
            timeline.addWidget(self.frame_spin)
            self.time_label = QtWidgets.QLabel(objectName="muted")
            self.time_label.setMinimumWidth(155)
            timeline.addWidget(self.time_label)
            self.loop_checkbox = QtWidgets.QCheckBox("Boucle")
            self.loop_checkbox.setToolTip("Reprendre au premier état à la fin de la lecture.")
            timeline.addWidget(self.loop_checkbox)
            timeline.addWidget(QtWidgets.QLabel("FPS cible", objectName="muted"))
            self.frame_rate = QtWidgets.QSpinBox()
            self.frame_rate.setRange(1, 60)
            self.frame_rate.setValue(60)
            self.frame_rate.setToolTip("Cadence visée. Le débit réel dépend du chargement et du rendu des données ; chaque état est affiché.")
            self.frame_rate.setFixedWidth(60)
            self.frame_rate.valueChanged.connect(self.change_frame_rate)
            timeline.addWidget(self.frame_rate)
            layout.addLayout(timeline)
            vector_row = QtWidgets.QHBoxLayout()
            vector_row.setSpacing(12)
            self.vector_checkbox = QtWidgets.QCheckBox("Vecteurs vitesse")
            self.vector_checkbox.setChecked(vectors)
            vector_row.addWidget(self.vector_checkbox)
            vector_row.addWidget(QtWidgets.QLabel("Flèches max.", objectName="muted"))
            self.vector_density = QtWidgets.QSpinBox()
            self.vector_density.setRange(50, 3000)
            self.vector_density.setSingleStep(100)
            self.vector_density.setValue(600)
            self.vector_density.setFixedWidth(85)
            vector_row.addWidget(self.vector_density)
            vector_row.addWidget(QtWidgets.QLabel("Échelle", objectName="muted"))
            self.vector_scale = QtWidgets.QDoubleSpinBox()
            self.vector_scale.setRange(.1, 4)
            self.vector_scale.setSingleStep(.25)
            self.vector_scale.setValue(1)
            self.vector_scale.setFixedWidth(85)
            self.vector_scale.setToolTip("Longueur proportionnelle à la vitesse ; échelle d'affichage réglable.")
            vector_row.addWidget(self.vector_scale)
            vector_row.addStretch()
            self.vector_info = QtWidgets.QLabel(objectName="muted")
            vector_row.addWidget(self.vector_info)
            layout.addLayout(vector_row)
            self.vector_checkbox.toggled.connect(self.refresh_vectors)
            self.vector_density.valueChanged.connect(self.refresh_vectors)
            self.vector_scale.valueChanged.connect(self.refresh_vectors)
            if headless:
                self.flow_widget = QtWidgets.QLabel()
                self.flow_widget.setMinimumSize(600, 380)
                self.flow_widget.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.plotter = pv.Plotter(off_screen=True, window_size=(1100, 570))
            else:
                from pyvistaqt import QtInteractor
                self.plotter = QtInteractor(flow, auto_update=False)
                self.flow_widget = self.plotter.interactor
            layout.addWidget(self.flow_widget, 1)
            self.probe_hint = "Survoler le champ pour lire les valeurs de la cellule pointée"
            self.probe_label = QtWidgets.QLabel(self.probe_hint, objectName="muted")
            layout.addWidget(self.probe_label)
            center.addWidget(flow, 1)
            curve = QtWidgets.QFrame(objectName="curveCard")
            layout = QtWidgets.QVBoxLayout(curve)
            layout.setContentsMargins(18, 12, 18, 6)
            row = QtWidgets.QHBoxLayout()
            self.curve_title = QtWidgets.QLabel("Convergence", objectName="title")
            row.addWidget(self.curve_title)
            row.addStretch()
            self.curve_summary = QtWidgets.QLabel(objectName="muted")
            row.addWidget(self.curve_summary)
            self.curve_combo = QtWidgets.QComboBox()
            self.curve_combo.addItems(["Convergence", "Distribution à la paroi", "Coefficients aérodynamiques"])
            self.curve_combo.currentIndexChanged.connect(self.change_curve)
            row.addWidget(self.curve_combo)
            curve_export = QtWidgets.QPushButton("Exporter")
            curve_export.setToolTip("Enregistrer le graphe affiché (PNG, SVG, PDF) ou ses données (CSV).")
            curve_export.clicked.connect(self.save_curve)
            row.addWidget(curve_export)
            layout.addLayout(row)
            self.figure = Figure(figsize=(10, 2), facecolor=THEME["card"], layout="constrained")
            self.canvas = FigureCanvasQTAgg(self.figure)
            self._history_background = None
            self.history_marker = None
            self.canvas.mpl_connect("draw_event", self.cache_history_background)
            self.canvas.mpl_connect("motion_notify_event", self.hover_curve)
            self.canvas.mpl_connect("figure_leave_event", self.hover_curve)
            self.curve_data = []  # (axes, x name, y name, x, y) of every curve on screen, for the hover readout.
            self.curve_highlights = {}
            self._curve_summary = ""
            self._coefficient_cache = (None, None)
            self.canvas.setFixedHeight(205)
            layout.addWidget(self.canvas)
            center.addWidget(curve)
            footer = QtWidgets.QLabel("  Coordonnées du maillage réel · Zoom à la molette · Déplacement à la souris · Espace : lecture/pause · F5 : actualiser", objectName="muted")
            footer.setContentsMargins(22, 0, 0, 10)
            outer.addWidget(footer)
            self.sequence = None
            self.frame_index = 0
            self.locked_ranges = {}
            self.play_timer = QtCore.QTimer(self)
            self.play_timer.setTimerType(QtCore.Qt.TimerType.PreciseTimer)
            self.play_timer.setSingleShot(True)
            self.play_timer.timeout.connect(self.advance_frame)
            self.playing = False
            self._next_frame_deadline = 0.0
            self.engine = None
            self._history_sequence = None
            self.theme_button.setText("Thème clair" if THEME["name"] == "sombre" else "Thème sombre")
            for keys, action in (("Space", self.toggle_playback), ("F5", self.reload), ("Ctrl+O", self.open_file), ("Ctrl+S", self.save_image)):
                QtGui.QShortcut(QtGui.QKeySequence(keys), self, activated=action)
            if self.plotter.iren is not None:
                self.plotter.iren.add_observer("EndInteractionEvent", self.refresh_vectors)
                self.plotter.iren.add_observer("MouseMoveEvent", self.hover_flow)
                self.plotter.iren.add_observer("LeaveEvent", lambda *args: self.probe_label.setText(self.probe_hint))
            self.set_solution(solution, field=field, history_override=history_path, requested_bounds=bounds, frame=frame)

        def set_solution(self, path, *, field="pressure", history_override=None, requested_bounds=None, frame=0, keep_view=False):
            sequence = SolutionSequence(path)
            if keep_view:  # Reload: the series may have grown or shrunk; a negative frame means the last one.
                frame = len(sequence) - 1 if frame < 0 else min(frame, len(sequence) - 1)
            mesh = sequence.read(frame)
            history_file = history_override or Path(path).with_name("residuals.csv")
            history = load_history(history_file) if history_file.exists() else None
            self.stop_playback()
            self.sequence = sequence
            self.history = history
            self.history_file = history_file
            duration = execution_time_text(history)
            self.execution_time_label.setText(duration)
            self.execution_time_caption.setText("Dernière mesure du calcul" if duration != "Indisponible" else "Durée non enregistrée dans l’historique")
            self.path = Path(path)
            if not keep_view:
                self.locked_ranges = {}
            self.file_label.setText(self.path.name)
            self.file_label.setToolTip(str(self.path.resolve()))
            self.frame_slider.blockSignals(True)
            self.frame_slider.setRange(0, len(sequence) - 1)
            self.frame_slider.setValue(frame)
            self.frame_slider.blockSignals(False)
            self.frame_spin.blockSignals(True)
            self.frame_spin.setRange(1, len(sequence))
            self.frame_spin.setValue(frame + 1)
            self.frame_spin.blockSignals(False)
            self.frame_slider.setEnabled(len(sequence) > 1)
            self.frame_spin.setEnabled(len(sequence) > 1)
            self.play_button.setEnabled(len(sequence) > 1)
            self.gif_button.setEnabled(len(sequence) > 1)
            self.display_frame(frame, field=field, bounds=requested_bounds, mesh=mesh, preserve_camera=keep_view and self.engine is not None)

        def apply_theme(self, name, save=False):
            """Switch the whole window, the flow view and the curves to another theme."""
            configure_theme(QtWidgets.QApplication.instance(), name)
            self.setStyleSheet(build_stylesheet())
            self.theme_button.setText("Thème clair" if name == "sombre" else "Thème sombre")
            self.figure.set_facecolor(THEME["card"])
            if self.engine is not None:
                self.engine.apply_theme()
                self._history_sequence = None  # The curves are redrawn with the new colours.
                self.draw_curve()
                if headless:
                    self.refresh_image()
            if save:
                QtCore.QSettings("CMFSN", "visualize").setValue("theme", name)

        def reload(self):
            """Read the same solution again, keeping the field, the camera and the locked colour ranges."""
            following_end = len(self.sequence) > 1 and self.frame_index == len(self.sequence) - 1
            try:
                self.set_solution(self.path, field=self.engine.current_field, history_override=self.history_file, frame=-1 if following_end else self.frame_index, keep_view=True)
            except (OSError, ValueError, RuntimeError) as error:
                QtWidgets.QMessageBox.critical(self, "Actualisation impossible", str(error))

        def display_frame(self, index, *, field=None, bounds=None, mesh=None, preserve_camera=True):
            mesh = self.sequence.read(index) if mesh is None else mesh
            field = field or (self.engine.current_field if self.engine is not None else "pressure")
            suppressed = self.plotter.suppress_rendering
            self.plotter.suppress_rendering = True
            try:
                if preserve_camera and self.engine is not None:
                    self.engine.set_mesh(mesh, field=field)
                else:
                    self.engine = SolutionViewer(mesh, field=field, bounds=bounds, plotter=self.plotter)
                self.frame_index = index
                self.frame_slider.blockSignals(True)
                self.frame_slider.setValue(index)
                self.frame_slider.blockSignals(False)
                self.frame_spin.blockSignals(True)
                self.frame_spin.setValue(index + 1)
                self.frame_spin.blockSignals(False)
                self.field_combo.blockSignals(True)
                if [self.field_combo.itemData(i) for i in range(self.field_combo.count())] != self.engine.fields:
                    self.field_combo.clear()
                    for name in self.engine.fields:
                        self.field_combo.addItem(FIELDS[name], name)
                self.field_combo.setCurrentIndex(self.engine.fields.index(field))
                self.field_combo.blockSignals(False)
                self.engine.set_colormap(self.palette_combo.currentText())
                self.engine.set_edges(self.edges.isChecked())
                self.apply_color_range()
                iteration = int(mesh.field_data["iteration"][0]) if "iteration" in mesh.field_data else None
                conditions = flow_summary(mesh)
                self.mesh_label.setText((f"{conditions}\n" if conditions else "") + f"{mesh.dimensions[0]-1:,} × {mesh.dimensions[1]-1:,} cellules\n{mesh.n_cells:,} cellules physiques" + (f"\nItération {iteration:,}" if iteration is not None else ""))
                time_value = self.sequence.times[index]
                self.time_label.setText(f"t = {time_value:.8g}\n{index + 1} / {len(self.sequence)} états" if np.isfinite(time_value) else "Solution unique")
                self.time_label.setToolTip("Temps physique dans les unités du calcul.")
                self.update_coefficients(mesh)
                self.draw_curve()
                self.refresh_vectors()
                self.update_stats()
            finally:
                self.plotter.suppress_rendering = suppressed
            if not suppressed:
                self.plotter.render()
                if headless:
                    self.refresh_image()

        def select_frame(self, index):
            try:
                self.display_frame(index)
            except (OSError, ValueError, RuntimeError) as error:
                self.stop_playback()
                self.frame_slider.blockSignals(True)
                self.frame_slider.setValue(self.frame_index)
                self.frame_slider.blockSignals(False)
                QtWidgets.QMessageBox.critical(self, "Lecture temporelle impossible", str(error))

        def stop_playback(self):
            self.playing = False
            self.play_timer.stop()
            self.play_button.setText("Lire")

        def schedule_next_frame(self):
            if self.playing:
                now = perf_counter()
                self._next_frame_deadline = max(self._next_frame_deadline + 1.0 / self.frame_rate.value(), now)
                self.play_timer.start(max(1, ceil((self._next_frame_deadline - now) * 1000)))

        def change_frame_rate(self, value):
            if self.playing:
                self._next_frame_deadline = perf_counter()
                self.schedule_next_frame()

        def toggle_playback(self):
            if self.playing:
                self.stop_playback()
            elif self.sequence is not None and len(self.sequence) > 1:
                if self.frame_index == len(self.sequence) - 1:
                    self.frame_slider.setValue(0)
                self.playing = True
                self.play_button.setText("Pause")
                self._next_frame_deadline = perf_counter()
                self.schedule_next_frame()

        def advance_frame(self):
            if self.frame_index + 1 < len(self.sequence):
                self.frame_slider.setValue(self.frame_index + 1)
                self.schedule_next_frame()
            elif self.loop_checkbox.isChecked():
                self.frame_slider.setValue(0)
                self.schedule_next_frame()
            else:
                self.stop_playback()

        def change_color_lock(self, checked):
            self.locked_ranges = {}
            if self.engine is not None:
                self.engine.set_field(self.engine.current_field)
                self.apply_color_range()
                self.update_stats()

        def apply_color_range(self):
            if self.color_lock.isChecked():
                limits = self.locked_ranges.setdefault(self.engine.current_field, self.engine.actor.mapper.scalar_range)
                self.engine.actor.mapper.scalar_range = limits
            self.plotter.render()

        def update_coefficients(self, mesh):
            result = aerodynamic_coefficients(mesh)
            if result is None:
                self.coefficients_label.setText("")
                self.coefficients_caption.setText("Coefficients aérodynamiques indisponibles")
                return
            self.coefficients_label.setText(f"Cl = {result['cl']:.5f}   ·   Cd = {result['cd']:.5f}   ·   Cm = {result['cm']:.5f}")
            self.coefficients_caption.setText(f"Cd : pression {result['cd_pressure']:.5f} + frottement {result['cd_friction']:.5f}" if result["viscous"] else "pression seule")

        def update_stats(self):
            self.field_title.setText(FIELDS[self.engine.current_field])
            self.min_label.setText(f"{self.engine.stats['min']:.6g}")
            self.max_label.setText(f"{self.engine.stats['max']:.6g}")
            invalid = self.engine.stats["invalid"]
            self.validity.setText("Toutes les valeurs sont finies" if invalid == 0 else f"{invalid:,} valeurs non finies")
            if headless:
                self.refresh_image()

        def change_field(self, index):
            if self.engine is not None and index >= 0:
                self.engine.set_field(self.field_combo.itemData(index))
                self.apply_color_range()
                self.update_stats()

        def change_palette(self, name):
            if self.engine is not None:
                self.engine.set_colormap(name)
                if headless:
                    self.refresh_image()

        def change_edges(self, visible):
            self.engine.set_edges(visible)
            if headless:
                self.refresh_image()

        def fit_domain(self):
            self.engine.reset_camera()
            self.refresh_vectors()
            if headless:
                self.refresh_image()

        def focus_airfoil(self):
            self.engine.focus((-.25, 1.25, -.65, .65))
            self.refresh_vectors()
            if headless:
                self.refresh_image()

        def refresh_vectors(self, *args):
            if self.engine is not None:
                self.engine.set_velocity_vectors(self.vector_checkbox.isChecked(), maximum=self.vector_density.value(), scale=self.vector_scale.value())
                active = self.vector_checkbox.isChecked()
                self.vector_density.setEnabled(active)
                self.vector_scale.setEnabled(active)
                self.vector_info.setText(f"{self.engine.vector_count} flèches visibles" if active else "")
                if headless:
                    self.refresh_image()

        def change_curve(self, index):
            self.curve_title.setText(self.curve_combo.currentText())
            self._history_sequence = None  # Either curve is drawn from scratch after a switch.
            if self.engine is not None:
                self.draw_curve()

        def draw_curve(self):
            (self.draw_history, self.draw_wall, self.draw_coefficients)[self.curve_combo.currentIndex()]()

        def set_curve_summary(self, text):
            self._curve_summary = text
            self.curve_summary.setText(text)

        def reset_curves(self):
            """Start a figure from scratch: no cached background, no hover target left from the previous one."""
            self.history_marker = None
            self._history_background = None
            self.curve_data = []
            self.curve_highlights = {}
            self.figure.clear()

        def hover_curve(self, event):
            """Show the values of the plotted point nearest to the cursor, and circle it."""
            target = None
            for axes, x_name, y_name, x, y in self.curve_data:
                if getattr(event, "inaxes", None) is axes and len(x):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        pixels = axes.transData.transform(np.column_stack([x, y]))
                    distance = np.hypot(pixels[:, 0] - event.x, pixels[:, 1] - event.y)
                    if np.isfinite(distance).any():
                        index = int(np.nanargmin(distance))
                        target = (axes, f"{x_name} = {x[index]:.6g}  ·  {y_name} = {y[index]:.6g}", x[index], y[index])
            for axes, highlight in self.curve_highlights.items():
                highlight.set_visible(target is not None and axes is target[0])
            if target is None:
                self.curve_summary.setText(self._curve_summary)
            else:
                axes, text, x, y = target
                if axes not in self.curve_highlights:
                    limits = axes.get_xlim(), axes.get_ylim()
                    self.curve_highlights[axes], = axes.plot([x], [y], linestyle="none", marker="o", markersize=9, markerfacecolor="none", markeredgecolor=THEME["marker"], markeredgewidth=1.6)
                    axes.set_xlim(limits[0]); axes.set_ylim(limits[1])
                self.curve_highlights[axes].set_data([x], [y])
                self.curve_summary.setText(text)
            self.canvas.draw_idle()

        def hover_flow(self, *args):
            """Read the cell under the cursor in the flow view."""
            if self.engine is None:
                return
            x, y = self.plotter.iren.get_event_position()
            renderer = self.plotter.renderer
            renderer.SetDisplayPoint(x, y, 0)
            renderer.DisplayToWorld()
            wx, wy, _, w = renderer.GetWorldPoint()
            self.show_probe(wx / w, wy / w) if w else self.probe_label.setText(self.probe_hint)

        def show_probe(self, x, y):
            self.probe_label.setText(self.engine.probe(x, y) or "Hors du maillage")

        def coefficient_history(self):
            """(x name, x, {cl, cd, cm}, source) from the solver history, else from the states of the series; None if neither."""
            names = [name for name, _ in COEFFICIENTS]
            if self.history is not None and set(names) <= set(self.history.dtype.names):
                use_time = "physical_time" in self.history.dtype.names and np.ptp(self.history["physical_time"]) > 0
                return ("t" if use_time else "Itération"), self.history["physical_time" if use_time else "iteration"], {name: self.history[name] for name in names}, "historique du calcul"
            if len(self.sequence) < 2:
                return None
            if self._coefficient_cache[0] is not self.sequence:
                indices = np.unique(np.linspace(0, len(self.sequence) - 1, min(len(self.sequence), MAX_COEFFICIENT_SAMPLES)).round().astype(int))
                values = {name: np.full(indices.size, np.nan) for name in names}
                reader = SolutionSequence(self.sequence.path)  # Separate reader: the displayed state is left alone.
                QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
                try:
                    for position, index in enumerate(indices):
                        result = aerodynamic_coefficients(reader.read(int(index)))
                        if result is not None:
                            for name in names:
                                values[name][position] = result[name]
                finally:
                    QtWidgets.QApplication.restoreOverrideCursor()
                self._coefficient_cache = (self.sequence, (self.sequence.times[indices], values))
            times, values = self._coefficient_cache[1]
            return "t", times, values, "états de la série" + (f" (1 sur {len(self.sequence) / times.size:.3g})" if times.size < len(self.sequence) else "")

        def draw_coefficients(self):
            """Evolution of Cl, Cd and Cm, each on its own scale."""
            self.reset_curves()
            data = self.coefficient_history()
            if data is None:
                axes = self.figure.add_subplot(111)
                axes.text(.5, .5, "Évolution indisponible : il faut une série temporelle (.pvd)\nou les colonnes cl, cd, cm dans residuals.csv", ha="center", va="center", transform=axes.transAxes, color=THEME["muted"])
                self.style_axes(axes, "", "")
                self.set_curve_summary("")
            else:
                x_name, x, values, source = data
                current = self.sequence.times[self.frame_index]
                for position, (name, title) in enumerate(COEFFICIENTS):
                    axes = self.figure.add_subplot(1, 3, position + 1)
                    axes.plot(x, values[name], color=THEME["curve"], linewidth=1.6, marker="o", markersize=3)
                    if x_name == "t" and np.isfinite(current):
                        axes.axvline(current, color=THEME["marker"], linewidth=1, linestyle="--")
                    self.style_axes(axes, "Temps physique" if x_name == "t" else "Itération", title)
                    finite = np.asarray(values[name], dtype=float)[np.isfinite(values[name])]
                    if finite.size and np.ptp(finite) < 1e-4:  # Round-off only: show it as the flat line it is.
                        centre = float(np.mean(finite))
                        axes.set_ylim(centre - 5e-4, centre + 5e-4)
                    axes.ticklabel_format(axis="y", useOffset=False)
                    self.curve_data.append((axes, x_name, title, np.asarray(x, dtype=float), np.asarray(values[name], dtype=float)))
                self.set_curve_summary(f"{len(x)} points · {source}")
            self.canvas.draw()

        def style_axes(self, axes, xlabel, ylabel):
            axes.set_facecolor(THEME["card"])
            axes.set_xlabel(xlabel, color=THEME["muted"], fontsize=9)
            axes.set_ylabel(ylabel, color=THEME["muted"], fontsize=9)
            axes.grid(True, which="major", color=THEME["grid"], linewidth=.7)
            axes.tick_params(colors=THEME["muted"], labelsize=8)
            for spine in axes.spines.values():
                spine.set_color(THEME["spine"])

        def draw_wall(self):
            """Cp against x in the cells touching the wall, for the displayed state."""
            self.reset_curves()
            name, x, values = wall_distribution(self.engine.mesh)
            axes = self.figure.add_subplot(111)
            axes.plot(x, values, color=THEME["curve"], linewidth=1.2, marker="o", markersize=2.5)
            if name == "cp":
                axes.invert_yaxis()  # Aerodynamic convention: suction upwards.
            self.style_axes(axes, "x", "Cp" if name == "cp" else "Pression")
            self.curve_data.append((axes, "x", "Cp" if name == "cp" else "p", np.asarray(x, dtype=float), np.asarray(values, dtype=float)))
            self.set_curve_summary(f"{x.size} cellules · première rangée à la paroi" + (" · axe Cp inversé" if name == "cp" else ""))
            self.canvas.draw()

        def draw_history(self):
            if getattr(self, "_history_sequence", None) is self.sequence:
                if self.history_marker is not None:
                    time_value = self.sequence.times[self.frame_index]
                    self.history_marker.set_xdata([time_value, time_value])
                    axes = self.history_marker.axes
                    if self._history_background is not None and tuple(axes.bbox.bounds) == self._history_background_bounds:
                        self.canvas.restore_region(self._history_background)
                        axes.draw_artist(self.history_marker)
                        self.canvas.blit(axes.bbox)
                    else:
                        self.canvas.draw_idle()
                return
            self._history_sequence = self.sequence
            self.reset_curves()
            axes = self.figure.add_subplot(111)
            axes.callbacks.connect("xlim_changed", self.invalidate_history_background)
            axes.callbacks.connect("ylim_changed", self.invalidate_history_background)
            if self.history is not None:
                use_time = "physical_time" in self.history.dtype.names and np.ptp(self.history["physical_time"]) > 0
                horizontal = self.history["physical_time"] if use_time else self.history["iteration"]
                axes.plot(horizontal, self.history["rms"], color=THEME["curve"], linewidth=1.8, marker="o", markersize=3.5)
                self.curve_data.append((axes, "t" if use_time else "Itération", "RMS", np.asarray(horizontal, dtype=float), np.asarray(self.history["rms"], dtype=float)))
                if use_time and self.sequence is not None and np.isfinite(self.sequence.times[self.frame_index]):
                    self.history_marker = axes.axvline(self.sequence.times[self.frame_index], color=THEME["marker"], linewidth=1, linestyle="--", animated=self.canvas.supports_blit)
                logarithmic = bool(np.all(self.history["rms"] > 0))
                if logarithmic:
                    axes.set_yscale("log")
                self.set_curve_summary(f"{self.history.size} mesures · " + ("échelle logarithmique" if logarithmic else "échelle linéaire"))
            else:
                axes.text(.5, .5, "Historique de convergence non fourni", ha="center", va="center", transform=axes.transAxes, color=THEME["muted"])
                self.set_curve_summary("")
            self.style_axes(axes, "Temps physique (unités du calcul)" if self.history is not None and use_time else "Itération", "RMS spatial")
            self.canvas.draw()

        def invalidate_history_background(self, *args):
            self._history_background = None

        def cache_history_background(self, event):
            if event.canvas is not self.canvas:
                return  # A file export (SVG, PDF) draws on its own canvas: nothing to cache for the screen.
            self._history_background = None
            if self.history_marker is not None and self.history_marker.get_animated():
                axes = self.history_marker.axes
                self._history_background = self.canvas.copy_from_bbox(axes.bbox)
                self._history_background_bounds = tuple(axes.bbox.bounds)
                axes.draw_artist(self.history_marker)

        def rendered_image(self):
            image = np.ascontiguousarray(self.plotter.screenshot(return_img=True))
            height, width, channels = image.shape
            fmt = QtGui.QImage.Format.Format_RGBA8888 if channels == 4 else QtGui.QImage.Format.Format_RGB888
            return QtGui.QImage(image.data, width, height, image.strides[0], fmt).copy()

        def refresh_image(self):
            if self.plotter.suppress_rendering:
                return
            image = self.rendered_image()
            self.flow_widget.setPixmap(QtGui.QPixmap.fromImage(image).scaled(self.flow_widget.size(), QtCore.Qt.AspectRatioMode.KeepAspectRatio, QtCore.Qt.TransformationMode.SmoothTransformation))

        def capture(self, output):
            QtWidgets.QApplication.processEvents()
            pixmap = self.grab()
            painter = QtGui.QPainter(pixmap)
            rect = QtCore.QRect(self.flow_widget.mapTo(self, QtCore.QPoint(0, 0)), self.flow_widget.size())
            painter.fillRect(rect, QtGui.QColor(THEME["plot"]))
            image = self.rendered_image().scaled(rect.size(), QtCore.Qt.AspectRatioMode.KeepAspectRatio, QtCore.Qt.TransformationMode.SmoothTransformation)
            painter.drawImage(rect.x()+(rect.width()-image.width())//2, rect.y()+(rect.height()-image.height())//2, image)
            painter.end()
            output.parent.mkdir(parents=True, exist_ok=True)
            if not pixmap.save(str(output)):
                raise OSError(f"Impossible d'écrire {output}")

        def export_curve(self, output):
            """The displayed graph as an image (from the extension), or its plotted points as CSV."""
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            if output.suffix.lower() == ".csv":
                if not self.curve_data:
                    raise ValueError("Le graphe affiché ne contient aucune donnée")
                with output.open("w", encoding="utf-8", newline="") as handle:
                    handle.write("courbe,x_nom,x,y\n")
                    for _, x_name, y_name, xs, ys in self.curve_data:
                        for x, y in zip(xs, ys):
                            handle.write(f"{y_name},{x_name},{x:.17g},{y:.17g}\n")
                return
            for highlight in self.curve_highlights.values():
                highlight.set_visible(False)  # The hover circle is not part of the graph.
            self.curve_summary.setText(self._curve_summary)
            import matplotlib
            # SVG keeps the labels as real text: they stay sharp at any zoom and can be read or searched in the file.
            with matplotlib.rc_context({"svg.fonttype": "none"}):
                self.figure.savefig(output, dpi=200, facecolor=self.figure.get_facecolor())
            self.canvas.draw_idle()

        def save_curve(self):
            name = CURVES[self.curve_combo.currentIndex()]
            formats = ["Image PNG (*.png)", "Image vectorielle SVG (*.svg)", "Document PDF (*.pdf)", "Données CSV (*.csv)"]
            # The wall Cp is exported as SVG by default: zooming into a vector image keeps every Cp value and tick readable.
            preferred = formats[1] if name == "paroi" else formats[0]
            suffix = ".svg" if name == "paroi" else ".png"
            chosen, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Exporter le graphe", str(self.path.with_name(f"{name}{suffix}")), ";;".join(formats), preferred)
            if chosen:
                try:
                    self.export_curve(Path(chosen))
                except (OSError, ValueError, RuntimeError) as error:
                    QtWidgets.QMessageBox.critical(self, "Export impossible", str(error))

        def export_gif(self, output, *, fps=10, stride=1, progress=None):
            """Animate the flow view over the series with a common colour range; the displayed state is restored afterwards."""
            if self.sequence is None or len(self.sequence) < 2:
                raise ValueError("Le GIF demande une série temporelle (.pvd) d'au moins deux états")
            self.stop_playback()
            field, start = self.engine.current_field, self.frame_index
            frames = gif_frames(len(self.sequence), stride=stride)
            reader = SolutionSequence(self.sequence.path)
            steps = 2 * frames.size
            low, high = np.inf, -np.inf
            for position, index in enumerate(frames):  # First pass: the range shared by every image.
                values = np.asarray(reader.read(int(index)).cell_data[field], dtype=float)
                values = values[np.isfinite(values)]
                if values.size:
                    low, high = min(low, float(values.min())), max(high, float(values.max()))
                if progress is not None and not progress(position + 1, steps):
                    return False
            if not low < high:
                low, high = low - 1e-6, high + 1e-6
            locked, ranges = self.color_lock.isChecked(), dict(self.locked_ranges)
            images = []
            self.color_lock.blockSignals(True)
            self.color_lock.setChecked(True)
            self.color_lock.blockSignals(False)
            self.locked_ranges = {field: (low, high)}
            try:
                for position, index in enumerate(frames):
                    self.display_frame(int(index))
                    time_value = self.sequence.times[index]
                    self.plotter.add_text(f"t = {time_value:.6g}" if np.isfinite(time_value) else f"État {index + 1}", name="gif_time", position="upper_left", color=THEME["text"], font_size=11)
                    self.plotter.render()
                    images.append(self.plotter.screenshot(return_img=True))
                    if progress is not None and not progress(frames.size + position + 1, steps):
                        return False
            finally:
                self.plotter.remove_actor("gif_time", render=False)
                self.color_lock.blockSignals(True)
                self.color_lock.setChecked(locked)
                self.color_lock.blockSignals(False)
                self.locked_ranges = ranges
                self.display_frame(start)
            write_gif(output, images, fps)
            return True

        def save_gif(self):
            chosen, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Exporter un GIF", str(self.path.with_name(f"{self.engine.current_field}.gif")), "Animation GIF (*.gif)")
            if not chosen:
                return
            dialog = QtWidgets.QProgressDialog("Préparation de l'animation…", "Annuler", 0, 100, self)
            dialog.setWindowTitle("Export du GIF")
            dialog.setMinimumDuration(0)
            dialog.setWindowModality(QtCore.Qt.WindowModality.WindowModal)

            def progress(done, total):
                dialog.setMaximum(total)
                dialog.setValue(done)
                dialog.setLabelText("Lecture des états…" if done <= total // 2 else "Rendu des images…")
                QtWidgets.QApplication.processEvents()
                return not dialog.wasCanceled()

            try:
                written = self.export_gif(Path(chosen), fps=self.frame_rate.value(), progress=progress)
            except (OSError, ValueError, RuntimeError) as error:
                QtWidgets.QMessageBox.critical(self, "Export impossible", str(error))
                return
            finally:
                dialog.close()
            if written:
                self.probe_label.setText(f"GIF enregistré : {chosen}")

        def save_image(self):
            chosen, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Exporter une image", str(self.path.with_suffix(".png")), "Image PNG (*.png)")
            if chosen:
                self.capture(Path(chosen))

        def open_file(self):
            chosen, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Ouvrir une solution", str(self.path.parent), "Solutions VTK (*.vts *.pvd)")
            if chosen:
                try:
                    self.set_solution(Path(chosen))
                except (OSError, ValueError, RuntimeError) as error:
                    QtWidgets.QMessageBox.critical(self, "Lecture impossible", str(error))

        def closeEvent(self, event):
            self.stop_playback()
            self.plotter.close()
            event.accept()

    return ViewerWindow()


def main(argv=None):
    import os
    parser = argparse.ArgumentParser(description="Visualiseur CMFSN : champs sur le maillage physique et convergence.")
    parser.add_argument("solution", nargs="?", type=Path)
    parser.add_argument("--field", choices=FIELDS, default="pressure")
    parser.add_argument("--history", type=Path)
    parser.add_argument("--frame", type=int, default=0, help="Indice temporel à afficher, à partir de zéro.")
    parser.add_argument("--vectors", action="store_true", help="Afficher le champ vectoriel de vitesse.")
    parser.add_argument("--bounds", type=float, nargs=4, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    parser.add_argument("--theme", choices=list(THEMES), help="Thème de l'interface (par défaut : le dernier choisi, sinon sombre).")
    parser.add_argument("--curve", choices=CURVES, default="convergence", help="Graphe affiché sous le champ.")
    parser.add_argument("--save-curve", type=Path, help="Enregistrer le graphe choisi (.png, .svg, .pdf) ou ses données (.csv), sans fenêtre.")
    parser.add_argument("--gif", type=Path, help="Enregistrer l'animation de la série (.pvd), sans fenêtre.")
    parser.add_argument("--fps", type=int, default=10, help="Cadence du GIF.")
    parser.add_argument("--stride", type=int, default=1, help="Un état sur N dans le GIF.")
    parser.add_argument("--screenshot", type=Path, help="Capture de l'interface sans ouvrir de fenêtre visible.")
    args = parser.parse_args(argv)
    os.environ.setdefault("QT_API", "pyside6")
    batch = bool(args.screenshot or args.save_curve or args.gif)
    if batch:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6 import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    configure_theme(app, args.theme or saved_theme())
    if args.solution is None:
        args.solution = latest_solution()
        if args.solution is not None:
            print(f"Solution la plus récente : {args.solution}")
    if args.solution is None:
        chosen, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Ouvrir une solution CMFSN", "", "Solutions VTK (*.vts *.pvd)")
        if not chosen:
            return 0
        args.solution = Path(chosen)
    window = create_window(args.solution, field=args.field, history_path=args.history, bounds=args.bounds, headless=batch, vectors=args.vectors, frame=args.frame)
    window.curve_combo.setCurrentIndex(CURVES.index(args.curve))
    window.show()
    app.processEvents()
    window.refresh_vectors()
    if batch:
        app.processEvents()
        if args.save_curve:
            window.export_curve(args.save_curve)
        if args.gif:
            window.export_gif(args.gif, fps=args.fps, stride=args.stride)
        if args.screenshot:
            window.capture(args.screenshot)
        window.close()
        return 0
    return app.exec()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Visualisation impossible : {error}", file=sys.stderr)
        raise SystemExit(1)
