# How to visualize protein motion from an MD trajectory

Use this workflow when you have an **aligned** protein structure and trajectory and want to identify flexible regions, inspect how their motion evolves over time, and make a structural visualization.

## Before you start

- Install [R](https://cran.r-project.org/) so the `Rscript` command is available. RMSX uses it to create plots and may install its plotting packages on the first run.
- Install ChimeraX only if you want to open a Flipbook visualization after analysis. VMD is an optional alternative renderer.
- Align the trajectory before analysis. RMSX measures internal residue motion; unremoved whole-protein translation or rotation will distort the result.

## Run the analysis

Install RMSX, then pass the topology followed by the trajectory:

```bash
pip install rmsx
rmsx topology.pdb trajectory.dcd --output_dir results --num_slices 9 --chain A --overwrite
```

Replace `topology.pdb`, `trajectory.dcd`, and chain `A` with your files and target chain. RMSX accepts topology and trajectory formats supported by MDAnalysis. Use `rmsx --help` for all options.

`--num_slices 9` divides the trajectory into nine time windows. Increase it for finer temporal detail when the trajectory contains enough frames. `--overwrite` replaces an existing managed output directory.

## Read the results

For a single selected chain, RMSX writes results to `results/chain_A_rmsx` (with `A` replaced by the selected chain). That directory contains:

- An RMSX CSV file and heatmap showing residue-level motion across time slices.
- `rmsd.csv` and `rmsf.csv` when RMSD/RMSF output is enabled.
- `slice_*_first_frame.pdb` structures. Their B-factor column holds the RMSX value for the corresponding slice, which allows structural coloring in molecular viewers.

RMSX identifies residues whose internal positional variation changes over the trajectory. High values point to regions with greater motion in that time slice; interpret them alongside the structure, simulation protocol, and conventional RMSD/RMSF outputs.

## Create a Flipbook visualization (optional)

To create and open a ChimeraX Flipbook from Python, use `run_rmsx_flipbook`:

```python
from rmsx import run_rmsx_flipbook

run_rmsx_flipbook(
    topology_file="topology.pdb",
    trajectory_file="trajectory.dcd",
    output_dir="results",
    num_slices=9,
    chain_sele="A",
    overwrite=True,
)
```

ChimeraX runs interactively, so a notebook cell may not finish until its window is closed. See the [method guide](choosing-a-method.md) for when a slice-wise RMSX map adds information beyond RMSF or RMSD.
