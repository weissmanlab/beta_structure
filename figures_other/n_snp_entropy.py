import json
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
n_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 25]
run_paths = ["runs/r001", "runs/r003", "runs/r008", "runs_structured/151"]  # interesting ones: 156, 234, 149
save_fig = True
###

fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C", "D", "D"],
        ["A", "A", "B", "B", "C", "C", "D", "D"],
    ],
    figsize = (8, 2),
    sharey = True,
    sharex = True
)

for label, run_path in zip(["A", "B", "C", "D"], run_paths):
    ax = axes[label]

    with open(run_path + "_entropy", "rb") as file:
        sample_entropy_all_n = pickle.load(file)
    with open(run_path + ".json", "r") as file:
        params = json.load(file)

    for sample, values in sample_entropy_all_n.items():
        sns.lineplot(x=n_vals, y=values, label=sample, ax=ax, legend=False)

    rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
    alpha = params["alpha"]
    if label == "D":
        ax.set_title(f"4B\nBeta ($\\alpha = {alpha:.3g}, \\rho = {rho:.3g}$)")
    else:
        ax.set_title(f"Kingman ($\\rho$ = {rho:.3g})")

    ax.set_xlim(0, 25)
    ax.set_ylim(0, 11)

    ax.text(-0.1, 1.1, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

fig.text(0.1, 0.5, "Entropy of n-SNP distribution (bits)", va="center", rotation="vertical")
fig.text(0.5, 0.00, "Number of strains n", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/manuscript/n_snp_entropy.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
