"""Figure style for the AMP MD chapter. Same locked page geometry and type scale as the
LARP6 chapter (fig_style.py), with AMP-specific condition colours (amp_config.CONDITIONS).

A4, 2.54 cm margins: text block 159.2 x 246.2 mm. Figures built at true size so they
place in Word at 100 percent. SVG keeps live text, PDF embeds TrueType: both open in
Illustrator as editable text.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

FULL_W_MM, HALF_W_MM = 159.2, 76.6
PT_PANEL, PT_NAME, PT_AXIS, PT_ANNOT, PT_FOOT = 11, 9, 7.5, 7, 6.5
C_GREY = "#BFBFBF"
C_SEG = ["#E8E8E8", "#F6F6F6"]          # alternating segment shading
# sequential (one hue, light to dark) and diverging (two hues, grey midpoint) ramps
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f4f8fd", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#a3421a", "#eb6834", "#f0c4ad", "#e6e6e6", "#b7d3f6", "#2a78d6", "#104281"])
HELIX_CMAP = ListedColormap(["#ededed", "#2a78d6"])


def apply():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": PT_AXIS, "axes.labelsize": PT_AXIS, "axes.titlesize": PT_NAME,
        "xtick.labelsize": PT_AXIS, "ytick.labelsize": PT_AXIS, "legend.fontsize": PT_ANNOT,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.0, "legend.frameon": False,
        "figure.dpi": 100, "savefig.dpi": 600, "savefig.facecolor": "white",
    })


def figure(width_mm=FULL_W_MM, height_mm=100.0):
    return plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))


def panel_letter(ax, letter, dx=None, dy=None):
    """Register a panel letter. It is placed at save time just outside the top-left corner of the
    panel's full extent (tick labels and axis label included), so it never overlaps the y label."""
    ax.figure.__dict__.setdefault("_letters", []).append((ax, letter))


def _place_letters(fig):
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, letter in getattr(fig, "_letters", []):
        bb = ax.get_tightbbox(r)
        x0, y1 = inv.transform((bb.x0, bb.y1))
        _, ytop = inv.transform((0, ax.get_window_extent(r).y1))
        fig.text(x0 - 0.004, max(y1, ytop), letter, fontsize=PT_PANEL, fontweight="bold", ha="right", va="top")


def residue_axis(ax, seq, segments, every=1, shade=True, label_segments=True):
    """x axis as one-letter sequence with residue numbers every 5, segments shaded."""
    n = len(seq)
    if shade:
        for i, (name, lo, hi) in enumerate(segments):
            ax.axvspan(lo - 0.5, hi + 0.5, color=C_SEG[i % 2], lw=0, zorder=0)
            if label_segments:
                ax.text((lo + hi) / 2, 1.0, name, transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=PT_FOOT, color="#555555")
    ax.set_xlim(0.5, n + 0.5)
    ax.set_xticks(range(1, n + 1, every))
    ax.set_xticklabels([f"{c}\n{i}" if (i % 5 == 0 or i == 1) else c for i, c in enumerate(seq, 1)][::every],
                       fontsize=PT_FOOT)
    ax.tick_params(axis="x", length=1.5, pad=1)


def save(fig, stem):
    _place_letters(fig)
    fig.savefig(stem + ".svg", format="svg", bbox_inches="tight")
    fig.savefig(stem + ".pdf", format="pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    return stem
