This repo contains simulation and analysis code for the manuscript:

Jiang, C. J., & Weissman, D. B. (2026). Bursts of reproduction can create genetic structure in frequently recombining bacterial populations. GENETICS. https://doi.org/10.1093/genetics/iyag132

# Project organization

Simulation parameters and seeds are stored as `.json` files.<br>

Scripts to generate figures in `figures_manuscript/`
- Fig 1: `example_tree.py`
- Fig 2: `kingman_progression.py`
  - Runs in `runs/*`
- Fig 3: `r_d_distribution_arity.py`
  - Runs simulated with `cluster/mass_sim_arity.py`
- Fig 4: `structure_platter.py`
  - Runs in `runs_structured/*`
- Fig 5: `three_peaks_segments.py`
  - Values in subplots D, E, and F computed using `pair_segments.py`
- Fig 6: `liu_and_good.py`
  - Runs in `runs/*` and `runs_structured/151.json`
  - Fraction of identical blocks computed using `cluster/frac_iden_blk.py`

Additionally, many scripts were reused across many figures for common tasks:
- Simulations run using `cluster/sim.py`.<br>
- Pairwise distances computed using `cluster/dist.py`
- $$\bar r_d$$ computed using `cluster/rd.py`
- Fraction of a pair's genome clonally inherited using `cluster/frac_clonal.py`

Some extra scripts for exploratory/miscellaneous plotting in  `./plt_*` and `figures_others/`<br>
