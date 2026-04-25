import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
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
paths = [
    "runs_structured/151_peak-1_pair_segments",
    "runs_structured/151_peak-2_pair_segments",
    "runs_structured/151_peak-3_pair_segments",
]
peaks = [
    "Fig. 4B left peak",
    "Fig. 4B middle peak",
    "Fig. 4B right peak",
]
order = ["Clonal",
         "Lineage 1\n recombined",
         "Lineage 2\n recombined",
         "Doubly\nrecombined"]
# palette (color unusued by other figures)
base_palette = sns.color_palette()[2:5]
base_palette.insert(2, "pink")
color_dict = {t: base_palette[i] for i, t in enumerate(order)}
###

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [["D", "D", "E", "E", "F", "F"],
     ["D", "D", "E", "E", "F", "F"]],
    figsize=(7, 2), 
    sharey=True,
)
plt.subplots_adjust(wspace=0.65)

bin_edges = np.linspace(0, 0.1, 30)

for label, peak, input_path in zip(["D", "E", "F"], peaks, paths):
    with open(input_path, "rb") as f:
        tmrcas = pickle.load(f)

    divergences = [t[1]*2*0.025 for t in tmrcas]
    recomb_types = []
    singly_clonal_types = {}
    for t in tmrcas:
        if t[0] == "Doubly":
            recomb_types.append("Doubly\nrecombined")
        elif t[0] == "Clonal":
            recomb_types.append("Clonal")
        else:
            if t[2] not in singly_clonal_types:
                singly_clonal_types[t[2]] = "Lineage 1\n recombined" if len(singly_clonal_types) == 0 else "Lineage 2\n recombined"
            recomb_types.append(f"{t[0]} {singly_clonal_types[t[2]]}".replace("Singly ", ""))

    ax = axes[label]
    sns.histplot(
        x=divergences, hue=recomb_types, hue_order=order,
        bins=bin_edges, multiple="stack", stat="probability",
        palette=color_dict, legend=(label=="D"), ax=ax
    )
    # for patch in ax.patches:
    #     if patch.get_x() > 0.01:
    #         patch.set_facecolor("grey")

    ax.text(-0.1, 1.1, rf"$\textbf{{{label}}}$", transform=ax.transAxes, fontweight="bold", va="top", ha="left")
    ax.set_title(peak)
    ax.set_ylim(0, 0.82)
    ax.set_xlim(-0.004, 0.1)
    ax.set_xticks(np.arange(0, 0.11, 0.03))

fig.text(0.5, -0.04, "Pairwise genetic distance ($d$) of segments", ha="center")
if save_fig:
    plt.savefig("../figures/manuscript/151_segments_tmrcas.png", dpi=500)
else:
    plt.show()
