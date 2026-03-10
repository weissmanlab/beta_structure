import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
input_path = "runs_mass/mass_sim_arity.csv"
###

df = pd.read_csv("runs_mass/mass_sim_arity.csv", header=None, index_col=None)

df.columns = ["arity", "T", "r_d"]

print(df)

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A"],
        ["A", "A"],
    ],
    figsize = (3, 3),
    sharey = True,
)

sns.scatterplot(x = df["r_d"], y = df["T"], hue = df["arity"], palette = "flare", ax = axes["A"], s=25)

plt.ylim(0, 0.7)
# plt.xlim(-0.001, 0.03)
plt.xscale("log")
plt.xlabel("$r_d$ of run")
plt.ylabel("T of burst")
plt.legend(title="Lineages captured\nin largest burst")

if save_fig:
    plt.savefig("../figures/manuscript/arity_vs_rd.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()

# plt.figure(figsize = (6,6))
# sns.scatterplot(x = df["r_d"], y = df["Lineages captured\nin largest burst"], s=55)
# plt.show()
