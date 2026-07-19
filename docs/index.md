# RMSX: visualize residue-level protein motion from MD trajectories

RMSX is a Python tool for finding, quantifying, and visualizing how regions of a protein move during a molecular-dynamics (MD) simulation. It turns an aligned topology and trajectory into residue-level RMSX maps, standard RMSD/RMSF summaries, and structures prepared for a ChimeraX Flipbook visualization.

## Use RMSX when

- You want to see which protein regions are flexible and how that pattern changes through a trajectory.
- You need publication-ready residue-by-time heatmaps plus structures colored with motion values.
- You want a concise ChimeraX visualization of an MD trajectory without manually preparing each time slice.

## Do not use RMSX when

- Your trajectory has not been aligned: remove overall translation and rotation first.
- You need to perform MD simulation, trajectory editing, or atom-level trajectory manipulation.
- You need collective-variable analysis or dimensionality reduction; consider PCA, tICA, or related workflows instead.

## Inputs and outputs

**Inputs:** an aligned topology file (for example PDB, PSF, or PRMTOP) and an MD trajectory supported by [MDAnalysis](https://www.mdanalysis.org/) (for example DCD, XTC, or TRR).

**Outputs:** RMSX heatmaps and CSV data, optional RMSD/RMSF tables and plots, PDB structures for each time slice with RMSX values in the B-factor field, and optional ChimeraX Flipbook files.

## Quick start

Install RMSX from PyPI:

```bash
pip install rmsx
```

RMSX uses **R** and an available `Rscript` command to make plots. ChimeraX or VMD is optional and only needed for Flipbook visualization.

```bash
rmsx topology.pdb trajectory.dcd --output_dir results --num_slices 9 --chain A --overwrite
```

See [Visualize protein motion](visualize-protein-motion.md) for a complete workflow and [Choose a method](choosing-a-method.md) to decide when RMSX is the right analysis.

## Citation

If you use RMSX in research, cite:

Beruldsen, F., de Freitas, M. V. & Antunes, D. A. *High resolution mapping of protein motions in time and space with RMSX and Flipbook.* **Scientific Reports** (2026). [https://doi.org/10.1038/s41598-026-39869-7](https://doi.org/10.1038/s41598-026-39869-7)
