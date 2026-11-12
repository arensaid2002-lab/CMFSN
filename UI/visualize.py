"""Inspect CMFSN structured-grid exports without interpolating cell values.

Install: python -m pip install -r UI/requirements.txt
Latest:  python UI/visualize.py                 (most recent results/run_*/solution)
Series:  python UI/visualize.py results/run_0001/solution.pvd
Single:  python UI/visualize.py results/run_0001/solution.vts
Zoom:    python UI/visualize.py results/run_0001/solution.vts --bounds -0.25 1.25 -0.75 0.75
Image:   python UI/visualize.py results/run_0001/solution.vts --screenshot pressure.png
Keys:    Space play/pause, F5 reload, Ctrl+O open, Ctrl+S export image
"""
from __future__ import annotations

import argparse
from functools import partial
from math import ceil
from time import perf_counter
from pathlib import Path
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
        self.plotter.set_background("white")
        self.actor = self.plotter.add_mesh(mesh, scalars=field, preference="cell", cmap="viridis", lighting=False, show_edges=False, nan_color="#cbd5e1", scalar_bar_args={"title": "", "vertical": False, "position_x": .22, "position_y": .025, "width": .56, "height": .10, "color": "#243248", "label_font_size": 12, "n_labels": 5})
        scalar_bar = self.plotter.scalar_bar
        scalar_bar.GetBackgroundProperty().SetColor(1.0, 1.0, 1.0)
        scalar_bar.SetDrawBackground(True)
        scalar_bar.GetLabelTextProperty().SetBackgroundColor(1.0, 1.0, 1.0)
        scalar_bar.GetLabelTextProperty().SetBackgroundOpacity(1.0)
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
        self.actor.prop.edge_color = "#66748a"
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


def configure_theme(app):
    """Use an explicit light palette instead of inheriting the desktop theme."""
    from PySide6 import QtCore, QtGui
    if app.property("cmfsnLightTheme"):
        return
    app.setStyle("Fusion")
    app.styleHints().setColorScheme(QtCore.Qt.ColorScheme.Light)
    palette = QtGui.QPalette()
    role = QtGui.QPalette.ColorRole
    colors = {
        role.Window: "#edf2f8", role.WindowText: "#243248", role.Base: "#ffffff",
        role.AlternateBase: "#f1f5f9", role.ToolTipBase: "#ffffff", role.ToolTipText: "#243248",
        role.Text: "#243248", role.Button: "#f8fafc", role.ButtonText: "#243248",
        role.BrightText: "#ffffff", role.Highlight: "#168b80", role.HighlightedText: "#ffffff",
        role.PlaceholderText: "#64748b", role.Accent: "#168b80",
    }
    for key, color in colors.items():
        palette.setColor(key, QtGui.QColor(color))
    for key in [role.WindowText, role.Text, role.ButtonText, role.PlaceholderText]:
        palette.setColor(QtGui.QPalette.ColorGroup.Disabled, key, QtGui.QColor("#64748b"))
    palette.setColor(QtGui.QPalette.ColorGroup.Disabled, role.Base, QtGui.QColor("#f1f5f9"))
    app.setPalette(palette)
    app.setProperty("cmfsnLightTheme", True)


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
            self.setStyleSheet("""
                QWidget { font-family: 'DejaVu Sans'; font-size: 13px; color: #243248; }
                QMainWindow, #root { background: #edf2f8; }
                #header { background: #14243c; }
                #brand { color: white; font-size: 24px; font-weight: 700; }
                #subtitle { color: #a6b8cf; font-size: 11px; }
                #filename { color: #cbd5e1; }
                #sidebar, #flowCard, #curveCard { background: white; border: 1px solid #e0e7f0; border-radius: 12px; }
                #section { color: #7b8ba1; font-size: 10px; font-weight: 700; }
                #title { font-size: 19px; font-weight: 600; color: #15253c; }
                #metric { font-size: 24px; font-weight: 600; color: #14243c; }
                #executionMetric { font-size: 20px; font-weight: 600; color: #14243c; }
                #muted { color: #75849a; font-size: 12px; }
                #badge { background: #e9f7f3; color: #18735c; border-radius: 6px; padding: 6px 10px; font-size: 10px; font-weight: 600; }
                QPushButton { border: 1px solid #dce5ef; border-radius: 7px; background: #f8fafc; padding: 9px 12px; font-weight: 600; }
                QPushButton:hover { background: #eef3f9; border-color: #a4b4c9; }
                QPushButton#primary { background: #168b80; color: white; border: none; }
                QPushButton#primary:hover { background: #11796f; }
                QSpinBox, QDoubleSpinBox { border: 1px solid #d9e2ef; border-radius: 5px; background: white; padding: 5px; }
                QComboBox { border: 1px solid #d9e2ef; border-radius: 7px; background: white; padding: 9px; }
                QComboBox::drop-down { width: 24px; border: none; }
                QCheckBox { spacing: 8px; padding: 6px 0; }
                QWidget:disabled { color: #64748b; }
                QPushButton:disabled { background: #f1f5f9; border-color: #dce5ef; color: #64748b; }
                QComboBox QAbstractItemView { background: white; color: #243248; selection-background-color: #168b80; selection-color: white; }
                QToolTip { background: white; color: #243248; border: 1px solid #cbd5e1; padding: 6px; }
                QSlider::groove:horizontal { background: #dce5ef; height: 6px; border-radius: 3px; }
                QSlider::sub-page:horizontal { background: #168b80; border-radius: 3px; }
                QSlider::handle:horizontal { background: #168b80; border: 1px solid #11796f; width: 13px; margin: -4px 0; border-radius: 6px; }
                QSlider::handle:horizontal:disabled { background: #94a3b8; border-color: #64748b; }
            """)
            root = QtWidgets.QWidget(objectName="root")
            outer = QtWidgets.QVBoxLayout(root)
            outer.setContentsMargins(0, 0, 0, 0)
            outer.setSpacing(0)
            self.setCentralWidget(root)
            header = QtWidgets.QFrame(objectName="header")
            top = QtWidgets.QHBoxLayout(header)
            top.setContentsMargins(24, 17, 24, 17)
            brand = QtWidgets.QVBoxLayout()
            brand.setSpacing(0)
            brand.addWidget(QtWidgets.QLabel("CMFSN", objectName="brand"))
            brand.addWidget(QtWidgets.QLabel("POST-TRAITEMENT  /  ÉCOULEMENT 2D", objectName="subtitle"))
            top.addLayout(brand)
            top.addStretch()
            self.file_label = QtWidgets.QLabel(objectName="filename")
            top.addWidget(self.file_label)
            reload_button = QtWidgets.QPushButton("Actualiser")
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
            self.palette_combo.addItems(["viridis", "plasma", "coolwarm", "cividis"])
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
            side.addWidget(QtWidgets.QLabel("Minimum", objectName="muted"))
            self.min_label = QtWidgets.QLabel(objectName="metric")
            side.addWidget(self.min_label)
            side.addWidget(QtWidgets.QLabel("Maximum", objectName="muted"))
            self.max_label = QtWidgets.QLabel(objectName="metric")
            side.addWidget(self.max_label)
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
            image_button.clicked.connect(self.save_image)
            side.addWidget(image_button)
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
            self.curve_combo.addItems(["Convergence", "Distribution à la paroi"])
            self.curve_combo.currentIndexChanged.connect(self.change_curve)
            row.addWidget(self.curve_combo)
            layout.addLayout(row)
            self.figure = Figure(figsize=(10, 2), facecolor="white", layout="constrained")
            self.canvas = FigureCanvasQTAgg(self.figure)
            self._history_background = None
            self.history_marker = None
            self.canvas.mpl_connect("draw_event", self.cache_history_background)
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
            for keys, action in (("Space", self.toggle_playback), ("F5", self.reload), ("Ctrl+O", self.open_file), ("Ctrl+S", self.save_image)):
                QtGui.QShortcut(QtGui.QKeySequence(keys), self, activated=action)
            if self.plotter.iren is not None:
                self.plotter.iren.add_observer("EndInteractionEvent", self.refresh_vectors)
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
            self.display_frame(frame, field=field, bounds=requested_bounds, mesh=mesh, preserve_camera=keep_view and self.engine is not None)

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
            if self.curve_combo.currentIndex() == 1:
                self.draw_wall()
            else:
                self.draw_history()

        def style_axes(self, axes, xlabel, ylabel):
            axes.set_facecolor("white")
            axes.set_xlabel(xlabel, color="#64748b", fontsize=9)
            axes.set_ylabel(ylabel, color="#64748b", fontsize=9)
            axes.grid(True, which="major", color="#e8edf4", linewidth=.7)
            axes.tick_params(colors="#64748b", labelsize=8)
            for spine in axes.spines.values():
                spine.set_color("#dce4ef")

        def draw_wall(self):
            """Cp against x in the cells touching the wall, for the displayed state."""
            self.history_marker = None
            self._history_background = None
            name, x, values = wall_distribution(self.engine.mesh)
            self.figure.clear()
            axes = self.figure.add_subplot(111)
            axes.plot(x, values, color="#168b80", linewidth=1.2, marker="o", markersize=2.5)
            if name == "cp":
                axes.invert_yaxis()  # Aerodynamic convention: suction upwards.
            self.style_axes(axes, "x", "Cp" if name == "cp" else "Pression")
            self.curve_summary.setText(f"{x.size} cellules · première rangée à la paroi" + (" · axe Cp inversé" if name == "cp" else ""))
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
            self.history_marker = None
            self._history_background = None
            self.figure.clear()
            axes = self.figure.add_subplot(111)
            axes.callbacks.connect("xlim_changed", self.invalidate_history_background)
            axes.callbacks.connect("ylim_changed", self.invalidate_history_background)
            if self.history is not None:
                use_time = "physical_time" in self.history.dtype.names and np.ptp(self.history["physical_time"]) > 0
                horizontal = self.history["physical_time"] if use_time else self.history["iteration"]
                axes.plot(horizontal, self.history["rms"], color="#168b80", linewidth=1.8, marker="o", markersize=3.5)
                if use_time and self.sequence is not None and np.isfinite(self.sequence.times[self.frame_index]):
                    self.history_marker = axes.axvline(self.sequence.times[self.frame_index], color="#d98b32", linewidth=1, linestyle="--", animated=self.canvas.supports_blit)
                logarithmic = bool(np.all(self.history["rms"] > 0))
                if logarithmic:
                    axes.set_yscale("log")
                self.curve_summary.setText(f"{self.history.size} mesures · " + ("échelle logarithmique" if logarithmic else "échelle linéaire"))
            else:
                axes.text(.5, .5, "Historique de convergence non fourni", ha="center", va="center", transform=axes.transAxes, color="#75849a")
                self.curve_summary.setText("")
            self.style_axes(axes, "Temps physique (unités du calcul)" if self.history is not None and use_time else "Itération", "RMS spatial")
            self.canvas.draw()

        def invalidate_history_background(self, *args):
            self._history_background = None

        def cache_history_background(self, event):
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
            painter.fillRect(rect, QtGui.QColor("white"))
            image = self.rendered_image().scaled(rect.size(), QtCore.Qt.AspectRatioMode.KeepAspectRatio, QtCore.Qt.TransformationMode.SmoothTransformation)
            painter.drawImage(rect.x()+(rect.width()-image.width())//2, rect.y()+(rect.height()-image.height())//2, image)
            painter.end()
            output.parent.mkdir(parents=True, exist_ok=True)
            if not pixmap.save(str(output)):
                raise OSError(f"Impossible d'écrire {output}")

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
    parser.add_argument("--screenshot", type=Path, help="Capture de l'interface sans ouvrir de fenêtre visible.")
    args = parser.parse_args(argv)
    os.environ.setdefault("QT_API", "pyside6")
    if args.screenshot:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6 import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    configure_theme(app)
    if args.solution is None:
        args.solution = latest_solution()
        if args.solution is not None:
            print(f"Solution la plus récente : {args.solution}")
    if args.solution is None:
        chosen, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Ouvrir une solution CMFSN", "", "Solutions VTK (*.vts *.pvd)")
        if not chosen:
            return 0
        args.solution = Path(chosen)
    window = create_window(args.solution, field=args.field, history_path=args.history, bounds=args.bounds, headless=args.screenshot is not None, vectors=args.vectors, frame=args.frame)
    window.show()
    app.processEvents()
    window.refresh_vectors()
    if args.screenshot:
        app.processEvents()
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
