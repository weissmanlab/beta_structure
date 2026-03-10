import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as cm
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
file_stem2 = "kingman_0.0"
file_stem3 = "kingman_1.5"
file_stem4 = "kingman_15"
file_stem5 = "kingman_45"
file_stem6 = "beta1.1_0.0"
file_stem7 = "beta1.1_1.5"
file_stem8 = "beta1.1_15"
file_stem9 = "beta1.1_45"
###

dir = "runs_mass/"
df2 = pd.read_csv(dir + file_stem2 + ".csv", header=None)
df3 = pd.read_csv(dir + file_stem3 + ".csv", header=None)
df4 = pd.read_csv(dir + file_stem4 + ".csv", header=None)
df5 = pd.read_csv(dir + file_stem5 + ".csv", header=None)
df6 = pd.read_csv(dir + file_stem6 + ".csv", header=None)
df7 = pd.read_csv(dir + file_stem7 + ".csv", header=None)
df8 = pd.read_csv(dir + file_stem8 + ".csv", header=None)
df9 = pd.read_csv(dir + file_stem9 + ".csv", header=None)

df2["model"] = r"Kingman"
df3["model"] = r"Kingman"
df4["model"] = r"Kingman"
df5["model"] = r"Kingman"
df6["model"] = r"Beta ($\alpha = 1.1$)"
df7["model"] = r"Beta ($\alpha = 1.1$)"
df8["model"] = r"Beta ($\alpha = 1.1$)"
df9["model"] = r"Beta ($\alpha = 1.1$)"

cdf1 = pd.concat([df2, df6], ignore_index=True)
cdf2 = pd.concat([df3, df7], ignore_index=True)
cdf3 = pd.concat([df4, df8], ignore_index=True)
cdf4 = pd.concat([df5, df9], ignore_index=True)
cdfs = [cdf1, cdf2, cdf3, cdf4]

titles = [r"$\rho = 0$", r"$\rho = 1.5$", r"$\rho = 15$", r"$\rho = 45$"]
indiv_runs_arrow = ["runs_structured/unstructured_beta", "runs_structured/151", "runs_structured/119"]

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A","A","A","A","B","B","B","B",".","E","E","E","E","E","E","E"],
        ["A","A","A","A","B","B","B","B",".","E","E","E","E","E","E","E"],
        ["A","A","A","A","B","B","B","B",".","E","E","E","E","E","E","E"],
        ["A","A","A","A","B","B","B","B",".","E","E","E","E","E","E","E"],
        ["C","C","C","C","D","D","D","D",".","E","E","E","E","E","E","E"],
        ["C","C","C","C","D","D","D","D",".","E","E","E","E","E","E","E"],
        ["C","C","C","C","D","D","D","D",".","E","E","E","E","E","E","E"],
        ["C","C","C","C","D","D","D","D",".","E","E","E","E","E","E","E"]
    ],
    figsize = (8, 4),
    sharex = False
)

fig.subplots_adjust(wspace=0.15, hspace=1.95, left = 0.15, bottom = 0.08)

for label, cdf, title, in zip(["A", "B", "C", "D"], cdfs, titles):

    ax = axes[label]

    x_col = cdf.columns[0]
    rcdf_data = []
    for model_name, group in cdf.groupby("model"):
        sorted_vals = np.sort(group[x_col])
        probs = 1-np.arange(1, len(sorted_vals)+1) / len(sorted_vals)
        rcdf_data.append(pd.DataFrame({
            x_col: sorted_vals,
            'reverse_cdf': probs,
            "model": model_name
        }))

    rcdf_df = pd.concat(rcdf_data, ignore_index=True)

    sns.scatterplot(
        data=rcdf_df, x=x_col, y="reverse_cdf", hue="model", ax=ax,
        legend=(label == "C"),
        edgecolor = "none",
        hue_order = [r"Kingman", r"Beta ($\alpha = 1.1$)"],
        palette = [sns.color_palette()[5], sns.color_palette()[3]],
        s = 7 
    )

    if label == "C":
        handles, labels = ax.get_legend_handles_labels()
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

        ax.legend(handles=new_handles, labels=labels, title="model")

    ax.text(-0.05, 1.1, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ### arrows to individual runs
    if label == "D":
        beta_df = rcdf_df[rcdf_df["model"] == r"Beta ($\alpha = 1.1$)"]

        for arrow_run, arrow_label in zip(indiv_runs_arrow, ["4A", "4B", "4C"]):
            with open(arrow_run + "_rd", "r") as file:
                arrow_r_d = float(file.read())

            x_target = arrow_r_d
            y_target = np.interp(x_target, beta_df[x_col], beta_df["reverse_cdf"])
            ax.annotate(
                arrow_label,
                xy=(x_target, y_target),
                xytext=(x_target + 0.12, y_target),
                arrowprops=dict(arrowstyle="->")
            )

    ax.set_xlim(0.0, 0.69)
    # ax.set_ylim(-0.01, 1.01)
    ax.set_yscale("log")
    # ax.set_xscale("log")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if ax != axes["A"] and ax != axes["C"]:
        ax.set_yticklabels([])
    if ax != axes["C"] and ax != axes["D"]:
        ax.set_xticklabels([])

### PLOT E
df = pd.read_csv("runs_mass/mass_sim_arity.csv", header=None, index_col=None)
df.columns = ["arity", "T", "r_d"]

base_cmap = mpl.colormaps["flare"]

resampling_indices = np.linspace(0, 1, 256)**(0.5) # sqrt transformation
new_flare = mcolors.ListedColormap(base_cmap(resampling_indices))

norm = mcolors.Normalize(
    vmin=df["r_d"].min(),
    vmax=df["r_d"].max(),
)

df = df.sort_values("r_d")

sns.scatterplot(
    x=df["arity"],
    y=df["T"],
    hue=df["r_d"],
    palette=new_flare,
    edgecolor="none",
    hue_norm=norm,
    ax=axes["E"],
    s=25,
    legend=False
)

sm = cm.ScalarMappable(norm=norm, cmap=new_flare)

cax = fig.add_axes([0.84, 0.58, 0.015, 0.25])
cbar = fig.colorbar(sm, cax=cax)

ticks=[0.0,0.1,0.2,0.3]
cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{t:.1f}" for t in ticks])

cbar.ax.set_title(r"$r_d$", rotation=0)

axes["E"].set_ylim(-0.01, 1.01)
# axes["E"].set_xscale("log")
axes["E"].set_xlabel("")
axes["E"].set_ylabel("Time of largest burst")
# axes["E"].legend(title="$\\bar r_d$")
axes["E"].set_title(r"$\rho = 45$, Beta ($\alpha = 1.1$)")

axes["E"].text(-0.05, 1.05, r"$\textbf{E}$", transform=axes["E"].transAxes, 
               fontweight="bold", va="top", ha="left")

### arrows to individual runs
for arrow_run, arrow_label in zip(indiv_runs_arrow, ["4A", "4B", "4C"]):
    with open(arrow_run + "_biggestburst", "r") as file:
        line = file.readline().strip()
        max_arity, biggest_burst_T = line.split(",")
        x_target = int(max_arity)
        y_target = float(biggest_burst_T)

    if arrow_label == "4A":
        axes["E"].annotate(
            arrow_label,
            xy=(x_target, y_target),
            xytext=(x_target + 3, y_target + 0.1),
            arrowprops=dict(arrowstyle="->")
        )
    elif arrow_label == "4B":
        axes["E"].annotate(
            arrow_label,
            xy=(x_target, y_target),
            xytext=(x_target + 3, y_target + 0.1),
            arrowprops=dict(arrowstyle="->")
        )
    else: 
        axes["E"].annotate(
            arrow_label,
            xy=(x_target, y_target),
            xytext=(x_target + 1, y_target + 0.1),
            arrowprops=dict(arrowstyle="->")
        )

sns.move_legend(axes["C"], "upper right")
fig.text(0.33, 0.00, "Normalized index of association ($\\bar r_d$)", ha="center")
fig.text(0.73, 0.00, "Size of largest burst (lineages captured)", ha="center")
fig.text(0.09, 0.5, "Probability", va="center", rotation="vertical")

if save_fig:
    plt.savefig("../figures/manuscript/r_d_distribution.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
