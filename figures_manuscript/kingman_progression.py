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
run_indices = ["r001", "r002", "r004", "r008"]
recomb_status_palette = {
    "Fully clonal": (1.0, 0.8784, 0.0),
    "Partially\nrecombined": sns.color_palette()[1],
    "Fully recombined": sns.color_palette()[0]
}
###

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C", "D", "D"],
        ["A", "A", "B", "B", "C", "C", "D", "D"],
    ],
    figsize = (8, 2),
    sharey = True,
)

def load_run(run_index):
    input_path = "runs/" + run_index
    with open(input_path + ".json", "r") as file:
        params = json.load(file)
    with open(input_path + "_dist", "rb") as file:
        dist = pickle.load(file)
    with open(input_path + "_rd", "r") as file:
        r_d = float(file.read())
    with open(input_path + "_frac_clonal", "rb") as file:
        clonal_tmrca = pickle.load(file)

    frac_clonal, clonal_tmrca = map(np.array, zip(*clonal_tmrca))
    clonal_tmrca = np.array([0 if x is None else x for x in clonal_tmrca])

    recomb_status = [
        "Fully recombined" if frac == 0 
        else "Partially\nrecombined" if 0 < frac < 1 
        else "Fully clonal" 
        for frac in frac_clonal
    ]

    return dist, recomb_status, r_d, params

# SUBPLOTS A-D
for label, run_index in zip(["A", "B", "C", "D"], run_indices):
    ax = axes[label]

    dist, recomb_status, r_d, params = load_run(run_index)

    sns.histplot(
        x=dist, stat="probability", hue=recomb_status,
        bins=40, multiple="stack", hue_order=["Fully clonal", "Partially\nrecombined", "Fully recombined"],
        palette=recomb_status_palette,
        ax=ax, legend = (label == "D")
    )

    ### inset MDS
    inset_ax = ax.inset_axes([0.05, 0.45, 0.35, 0.35])
    dist_matrix = squareform(dist)

    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_coords = mds.fit_transform(dist_matrix)

    sns.scatterplot(x=mds_coords[:, 0], y=mds_coords[:, 1],
                    ax=inset_ax, color=sns.color_palette()[5],
                    edgecolor="white", linewidth=0.2, s=12)

    inset_ax.set_ylim(-0.04, 0.04)
    inset_ax.set_ylim(-0.025, 0.025)

    inset_ax.set_xticklabels([])
    inset_ax.set_yticklabels([])
    inset_ax.set_xlabel("")
    inset_ax.set_ylabel("")
    ###

    # panel labels
    ax.text(-0.1, 1.1, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")
    
    # r_d value
    ax.text(0.05, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
            verticalalignment="top")
    
    # subplot title/labels
    rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
    ax.set_title(f"$\\rho = {rho:.4g}$")
    ax.set_xlabel("")
    if ax != axes["A"]:
        ax.set_ylabel("")
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

sns.move_legend(axes["D"], "lower left")

fig.text(0.5, 0.00, "Pairwise genetic distance ($d$)", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/manuscript/kingman_progression.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
