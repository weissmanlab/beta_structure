import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from sklearn.manifold import MDS
from scipy.spatial.distance import squareform
import scienceplots
plt.style.use("science")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.titlesize": 10,
})

###
save_fig = True
run_indices = ["unstructured_beta", "151", "119"]
recomb_status_palette = {
    "Partially\nrecombined": sns.color_palette()[1],
    "Fully\nrecombined": sns.color_palette()[0]
}
###

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "A", "B", "B", "C", "C"],
    ],
    figsize = (6, 2),
    sharex = True
)

def load_run(run_index):
    input_path = "runs_structured/" + run_index
    with open(input_path + ".json", "r") as file:
        params = json.load(file)
    with open(input_path + "_dist", "rb") as file:
        dist = pickle.load(file)
        dist = np.array(dist)
    with open(input_path + "_frac_clonal", "rb") as file:
        clonal_tmrca = pickle.load(file)
    with open(input_path + "_rd", "r") as file:
        r_d = float(file.read())

    frac_clonal, clonal_tmrca = map(np.array, zip(*clonal_tmrca))

    recomb_status = [
        "Fully\nrecombined" if frac == 0 
        else "Partially\nrecombined"
        for frac in frac_clonal
    ]

    return dist, recomb_status, r_d, params

bin_edges = np.linspace(0, 0.036, 40)

for label, run_index in zip(["A", "B", "C"], run_indices):
    dist, recomb_status, r_d, params = load_run(run_index)
    ax = axes[label]

    ### inset MDS & r_d value
    if label == "C":
        inset_ax = ax.inset_axes([0.60, 0.50, 0.35, 0.35])
        ax.text(0.50, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
                verticalalignment="top")
    else:
        inset_ax = ax.inset_axes([0.05, 0.50, 0.35, 0.35])
        ax.text(0.05, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
                verticalalignment="top")
    dist_matrix = squareform(dist)

    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_coords = mds.fit_transform(dist_matrix)

    jitter_strength = 0.0005
    jittered_coords = mds_coords + np.random.normal(loc=0, scale=jitter_strength, size=mds_coords.shape)

    sns.scatterplot(x=jittered_coords[:, 0], y=jittered_coords[:, 1],
                    ax=inset_ax, color=sns.color_palette()[5],
                    edgecolor="white", linewidth=0.2, s=12)

    inset_ax.set_ylim(-0.04, 0.04)
    inset_ax.set_ylim(-0.026, 0.026)

    inset_ax.set_xticklabels([])
    inset_ax.set_yticklabels([])
    inset_ax.set_xlabel("")
    inset_ax.set_ylabel("")

    ## main histogram
    sns.histplot(
        x=dist, stat="probability", hue=recomb_status,
        bins=bin_edges, multiple="stack", hue_order=["Partially\nrecombined", "Fully\nrecombined"],
        ax=ax, palette=recomb_status_palette, legend = (label == "A")
    )

    ax.text(-0.1, 1.10, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax.set_ylim(0, 0.75)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if label != "A":
        ax.set_yticklabels([])

    rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
    ax.set_title(f"$\\rho = {rho:.4g}$")

sns.move_legend(axes["A"], "lower left")
axes["A"].set_ylabel("Probability")
fig.text(0.5, 0.00, "Pairwise genetic distance ($d$)", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/manuscript/structure_platter.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
