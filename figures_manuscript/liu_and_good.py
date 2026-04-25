import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
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
kingman_indices = ["r001", "r002", "r003", "r004", "r008"]
beta_index = "151"
blk_size = 1000 # for analytical prediction
recomb_status_palette = {
    "Partially recombined": sns.color_palette()[1],
    "Fully recombined": sns.color_palette()[0]
}
###

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "A", "B", "B", "B", "C"],
        ["A", "A", "A", "B", "B", "B", "C"],
        ["A", "A", "A", "B", "B", "B", "C"]
    ],
    figsize = (7, 3),
    sharey = True,
)

## RECOMBINANT LINE
def expected_dist(f, avg_d):
    mu     = params["mu"]
    r      = params["r"]
    t      = params["tract_length"]

    # per base rate of replacement by recombination
    # R = r * (t) * np.exp(-blk_size/t)
    R = r * (t)

    denom = R + mu * blk_size

    term_recomb = avg_d * (1 - f**(R/denom)) # SNPs introduced by recombination
    term_mut    = f**(R/denom) * ((mu)*np.log(f))/denom # SNPs introduced by mutation

    return term_recomb - term_mut

plt.subplots_adjust(wspace=0.35)

# SUBPLOT E: KC
axes["A"].set_xlim(0, 1)
axes["A"].text(-0.05, 1.06, r"$\textbf{A}$", transform=axes["A"].transAxes, 
               fontweight="bold", va="top", ha="left")
axes["A"].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
fig.text(0.28, 0.91, "Kingman", ha="center")

palette = sns.color_palette("plasma", n_colors=5)
# palette = sns.color_palette("PuOr", n_colors=5)
# palette[2] = (160/255, 102/255, 88/255)
# palette[1] = palette[0]
# palette[0] = (220/255, 69/255, 0/255)
# palette[3] = palette[4]
# palette[4] = (205/255, 71/255, 120/255)
for kingman_index, color in zip(kingman_indices, palette):
    input_path = "runs/" + kingman_index
    with open(input_path + ".json", "r") as file:
        params = json.load(file)
    with open(input_path + "_frac_iden_blk", "rb") as file:
        frac_iden_blk = pickle.load(file)
    with open(input_path + "_dist", "rb") as file:
        dist = pickle.load(file)
    avg_d = np.mean(dist)

    # simulation points
    rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
    sns.scatterplot(x=frac_iden_blk, y=dist,
                    ax=axes["A"],
                    s=14, 
                    color = color,
                    label=f"$\\rho$ = {rho:.3g}")

    # expectation lines
    r_x = np.linspace(1e-10, 1, 1000)
    r_y = expected_dist(r_x, avg_d)
    axes["A"].plot(r_x, r_y, 
                   linestyle="dashed",
                   color="black", 
                   alpha=0.5)

# legend adjustments
sns.move_legend(axes["A"], "upper left")
handles, labels = axes["A"].get_legend_handles_labels()
new_handles = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label=label,
        markerfacecolor=handle.get_facecolor()[0],
        markersize=(28**0.5),
        linestyle='None'
    )
    for handle, label in zip(handles, labels)
]
axes["A"].legend(handles=new_handles, labels=labels)

### SUBPLOT A - B: BC
input_path = "runs_structured/" + beta_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
with open(input_path + "_frac_iden_blk", "rb") as file:
    frac_iden_blk = pickle.load(file)
with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)
    avg_d = np.mean(dist)
with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

frac_clonal, _ = map(np.array, zip(*clonal_tmrca))
recomb_status = [
    "Fully recombined" if frac == 0 
    else "Partially recombined"
    for frac in frac_clonal
]

# BC points
sns.scatterplot(x=frac_iden_blk, y=dist, hue=recomb_status, hue_order=["Partially recombined", "Fully recombined"], palette=recomb_status_palette, ax=axes["B"], s=14)

# legend adjustments 2
sns.move_legend(axes["B"], "upper left")
handles, labels = axes["B"].get_legend_handles_labels()
new_handles = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label=label,
        markerfacecolor=handle.get_markerfacecolor(),
        markersize=(28**0.5),
        linestyle='None'
    )
    for handle, label in zip(handles, labels)
]
axes["B"].legend(handles=new_handles, labels=labels)

new_handles.append(Line2D([0], [0], color='black', linestyle='dashed', alpha=0.5, label='Expected (recombinant)'))
labels.append('Expected (Kingman)')
new_handles.append(Line2D([0], [0], color='red', linestyle='dashed', alpha=0.5, label='Expected (burst-adjusted)'))
labels.append('Expected (Burst-adjusted)')

axes["B"].legend(handles=new_handles, labels=labels)

# expectation line
r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_dist(r_x, avg_d)
axes["B"].plot(r_x, r_y, 
               linestyle='dashed', 
               color="black",
               alpha=0.5
               )


S_burst, T_burst = open(input_path + "_largest_burst").read().strip().split(",")
S_burst = int(S_burst)
T_burst = float(T_burst)
def adjusted_expected_dist(f):
    mu     = params["mu"]
    avg_f = np.mean(frac_iden_blk)

    return avg_d - ((f - avg_f) * (avg_d - 2 * mu * T_burst))/(np.exp(-2*mu*blk_size*T_burst)-avg_f)

adjusted_r_x = np.linspace(1e-10, 1, 1000)
adjusted_r_y = adjusted_expected_dist(adjusted_r_x)
axes["B"].plot(adjusted_r_x, adjusted_r_y, 
               linestyle='dashed', 
               color="red",
               alpha=0.5
               )

# BC marginal histogram
sns.histplot(y = dist, bins=50, ax = axes["C"],hue=recomb_status, multiple="stack", hue_order=["Partially recombined", "Fully recombined"], palette=recomb_status_palette, stat="probability", legend=False)

# settings
axes["B"].set_yticks([0.0, 0.01, 0.02, 0.03, 0.04])
axes["B"].set_xlim(0, 1)
axes["B"].set_ylim(0, 0.042)
axes["B"].text(-0.05, 1.06, r"$\textbf{B}$", transform=axes["B"].transAxes, 
               fontweight="bold", va="top", ha="left")
axes["C"].set_xlabel("")
axes["C"].set_xticks([])

alpha = params["alpha"]
fig.text(0.63, 0.91, f"Run 4B -- Beta ($\\alpha = {alpha:.3g}$)", ha="center")
fig.text(0.5, -0.01, "Fraction of identical 1 kb blocks ($f$)", ha="center")
fig.text(0.06, 0.5, "Pairwise genetic distance ($d$)", rotation="vertical", va="center")

if save_fig:
    plt.savefig("../figures/manuscript/liu_and_good.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
