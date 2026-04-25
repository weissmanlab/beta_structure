import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "52"
###

input_path = "runs_structured/" + run_index
blk_size = 1000 # for analytical prediction

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_frac_iden_blk", "rb") as file:
    frac_iden_blk = pickle.load(file)

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)

S_burst, T_burst = open(input_path + "_largest_burst").read().strip().split(",")
S_burst = int(S_burst)
T_burst = float(T_burst)

avg_f = np.mean(frac_iden_blk)
avg_d = np.mean(dist)
print("Average pi:", avg_d)

### NULL POINTS
null_index = "r001"
null_path = "runs/" + null_index
with open(null_path + "_frac_iden_blk", "rb") as file:
    null_frac_iden_blk = pickle.load(file)
with open(null_path + "_dist", "rb") as file:
    null_dist = pickle.load(file)

## JOINT PLOT
g = sns.jointplot(
    x=frac_iden_blk, 
    y=dist, 
    height=6, 
    space=0,
    xlim=(0,1),
    ylim=(0,0.041),
    marginal_kws={"bins": 160}
)

## NULL POINTS
sns.scatterplot(x=null_frac_iden_blk, y=null_dist, color='grey')

## NULL LINE
n_x = np.linspace(1e-3, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
g.ax_joint.plot(n_x, n_y, color='grey', linestyle='--')

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

r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_dist(r_x, avg_d)
g.ax_joint.plot(r_x, r_y, color='red', linestyle='--')


def adjusted_expected_dist(f):
    mu     = params["mu"]

    # # ORIGINAl
    # avg_d_cond = 0.029
    # return avg_d_cond - (f*np.exp(2*mu*blk_size*T_burst))*(avg_d_cond - 2*mu*T_burst)
    # REANCHORED:
    return avg_d - ((f - avg_f) * (avg_d - 2 * mu * T_burst))/(np.exp(-2*mu*blk_size*T_burst)-avg_f)

adjusted_r_x = np.linspace(1e-10, 1, 1000)
adjusted_r_y = adjusted_expected_dist(adjusted_r_x)
g.ax_joint.plot(adjusted_r_x, adjusted_r_y, 
               linestyle='--', 
               color="black")

## labels
g.set_axis_labels("Proportion of 1kb base blocks identical", 
                  "Pairwise mean number of nucleotide differences", 
                  fontsize=12)
rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Fraction of identical blocks vs distance (" + model_str + ", $\\rho$=" + str(rho)  + ")")

if save_fig:
    g.figure.savefig("/Users/cjdjpj/Desktop/" + run_index + "_frac_iden_blk.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
    plt.show()
