import msprime
import scipy
import numpy as np

nsample = 10
l = 1

def draw_svg_func():
    return ts.draw_svg(size=(500, 300),
                       symbol_size = 0,
                       node_labels={},
                       x_axis=False,
                       style="""
        path.edge { stroke-width: 6px; stroke: black; }
        circle.node { display: none; }
    """
                       )

ts = msprime.sim_ancestry(
    samples=nsample,
    ploidy=1,
    sequence_length=l,
)

tree_svg = draw_svg_func()

with open("../figures/manuscript/kingman.svg", "w") as f:
    f.write(tree_svg)

def T2(a, N):
    """Returns the expected pairwise coalescence time in a Beta coalescent"""
    return np.power(1 + 1 / np.exp2(a - 1) / (a - 1), a) * np.power(N, a - 1) / a / scipy.special.beta(2 - a, a)

def n_beta(a, T2):
    """Returns the N necessary to make a Beta coalescent with the given alpha = a have the specified pairwise coalescence time"""
    """if returns 0, means N is unimportant. if returns inf, means T unattainable with given alpha"""
    return ((T2 * a * scipy.special.beta(2 - a, a)) / ((1 + 1 / (2**(a - 1) * (a - 1)))**a))**(1 / (a - 1))

ts = msprime.sim_ancestry(
    samples=nsample,
    model=msprime.BetaCoalescent(alpha = 1.1),
    population_size = n_beta(1.1, 1),
    ploidy=1,
    sequence_length=l,
)

tree_svg = draw_svg_func()

with open("../figures/manuscript/beta.svg", "w") as f:
    f.write(tree_svg)
