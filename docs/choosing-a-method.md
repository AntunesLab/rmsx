# RMSX, RMSF, RMSD, PCA, and normal modes: which should you use?

These approaches answer related but different questions. RMSX is most useful when you want a residue-level view of motion **over time**, not a single summary of the entire simulation.

| Method | Best question | Choose it when | RMSX is different because |
| --- | --- | --- | --- |
| **RMSX** | Which residues move, and how does that motion vary across trajectory windows? | You need a time-resolved residue map and structural visualization. | It summarizes residue-level variation separately for consecutive slices. |
| **RMSF** | Which residues are most flexible across the whole trajectory? | You want one conventional per-residue flexibility summary. | It can hide when a flexible state occurred. |
| **RMSD** | How far has the structure moved from a reference over time? | You want a global stability, equilibration, or drift check. | It is a global metric, not a residue-level motion map. |
| **PCA / tICA** | What collective conformational changes dominate the simulation? | You need low-dimensional states, transitions, or clustering. | It does not replace dimensionality-reduction analysis. |
| **Normal-mode analysis** | What motions are plausible near a reference structure? | You are studying predicted low-frequency modes without a trajectory. | It analyzes a model rather than measured MD frames. |

## A practical workflow

1. Align and quality-check the trajectory.
2. Inspect RMSD to understand equilibration and major drift.
3. Use RMSF for an overall flexibility summary.
4. Use RMSX when you need to locate the time windows in which specific residue regions become more or less mobile.
5. Use PCA/tICA or normal modes when the research question is about collective conformational states rather than a residue-level visualization.

RMSX complements these methods; it should not be presented as a replacement for every type of trajectory analysis. For a runnable analysis, see [Visualize protein motion](visualize-protein-motion.md).
