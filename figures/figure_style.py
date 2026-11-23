"""Common presentation settings for the archived-data figures."""

import matplotlib as mpl


GAP_COLOR = "red"
GAP_ALPHA = 0.85
ALLOWED_COLOR = "royalblue"
FLAT_COLOR = (0, 0.7, 0)
FIELD_CMAP = "viridis"
SHAPE_STYLE = {
    "triangle": ("Triangle", "#2a78d6", "^"),
    "square": ("Square", "#eda100", "o"),
    "pentagon": ("Pentagon", "#e87ba4", "p"),
    "hexagon": ("Hexagon", "#008300", "h"),
    "horizontal_rectangle": ("Horiz. rectangle", "#eb6834", "s"),
    "vertical_rectangle": ("Vert. rectangle", "#1baf7a", "D"),
}


def apply_style():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "legend.frameon": False,
        "grid.alpha": 0.3,
        "grid.linewidth": 0.5,
        "pdf.fonttype": 42,
    })
