# This note shows how to bolt a new analysis method onto the RMSX/Flipbook pipeline, using a per-residue lDDT time-series as the example.

Goals:

* Reuse the **existing plumbing** (RMSX, shift maps, Flipbook).
* Keep each new metric in a **small, self-contained module**.
* Produce the same style of outputs (**CSV + PDB B-factors**) so you get:

  * residue × slice heatmaps
  * Flipbook visualizations “for free”
* Internally compute **lDDT ∈ [0, 1]**, but **store and plot `1 − lDDT`**, so that:

  * **high values = more motion / less local structural agreement**,
  * consistent with RMSX and shift maps (“hotter = more unstable”).

---

## 1. RMSX vs Shift Map: Window-Based vs Template-Based

RMSX and shift maps both output a **ResidueID × slice** matrix, but the logic is fundamentally different:

* **RMSX = window-based**
  “How much does residue *i* fluctuate *within this time window*?”

* **Trajectory / shift map = template-based**
  “How far has residue *i* moved away from a *fixed reference structure* over time?”

### 1.1. How the pipelines differ in code

Think of a trajectory split into N slices:

```text
Frames:  [ 0 ...  99] [100 ... 199] [200 ... 299] ...
Slices:     slice 1      slice 2        slice 3
```

#### RMSX: window-based (RMSF per window)

* Uses **all frames in each slice**.
* For slice *k*, you compute RMSF over the window `[slice_start : slice_end)`.

Pattern (simplified from `run_rmsx()` → `process_trajectory_slices_by_num()`):

```python
def process_trajectory_slices_by_num(u, output_dir, total_size, num_slices, ...):
    # decide frame ranges for each slice
    for i in range(num_slices):
        slice_start = ...
        slice_end   = ...
        u.trajectory[slice_start]

        ca_atoms  = u.select_atoms("protein and name CA")
        rmsf_calc = RMSF(ca_atoms)
        rmsf_calc.run(start=slice_start, stop=slice_end)

        df_slice = pd.DataFrame(
            {f"slice_{i+1}.dcd": rmsf_calc.results.rmsf},
            index=[res.resid for res in ca_atoms.residues],
        )
        # concat df_slice into all_data
```

**Interpretation:**
For slice 2, the value at residue *i* is “how much residue *i* wiggles between frames 100–199”.

---

#### Shift map: template-based (distance to reference)

* Picks a **reference frame**, usually the first frame of the run.
* For slice *k*, you look at **one frame** (e.g. the first frame in the slice) and
  measure the distance to that reference.

Pattern (simplified from `run_shift_map()` → `process_trajectory_shifts_by_size()`):

```python
def process_trajectory_shifts_by_size(u, output_dir, total_size, slice_size, ...):
    # 1) Set reference frame
    u.trajectory[start_frame]
    atoms_full_ref = u.select_atoms(selection_str)
    atoms_calc_ref = atoms_full_ref.select_atoms("name CA")
    ref_coords     = atoms_calc_ref.positions.copy()

    # 2) For each slice, compare to reference
    adjusted_total_size = ...
    n_slices = adjusted_total_size // slice_size

    for i in range(n_slices):
        slice_index = start_frame + i * slice_size  # first frame in slice
        u.trajectory[slice_index]

        atoms_full_current = u.select_atoms(selection_str)
        atoms_calc_current = atoms_full_current.select_atoms("name CA")
        curr_coords        = atoms_calc_current.positions

        distances = np.linalg.norm(curr_coords - ref_coords, axis=1)

        df_slice = pd.DataFrame(
            {f"slice_{i+1}.dcd": distances},
            index=[res.resid for res in atoms_calc_current.residues],
        )
        # concat df_slice into all_data
```

**Interpretation:**
For slice 2, the value at residue *i* is “distance between residue *i* at frame 100 and residue *i* at frame 0”.

---

### 1.2. Concrete toy example

One residue, three slices, frames:

* Slice 1: 0–9
* Slice 2: 10–19
* Slice 3: 20–29

Assume the residue’s position behaves like this:

```text
Frames 0–9:    always at 0 Å
Frames 10–19:  oscillates between 0 and 2 Å (wiggly)
Frames 20–29:  always at 5 Å (shifted but rigid)
Reference frame for shift map: frame 0 (position = 0 Å)
```

**RMSX (window-based, RMSF):**

* Slice 1: RMSF ≈ 0 (no fluctuation)
* Slice 2: RMSF > 0 (wiggles)
* Slice 3: RMSF ≈ 0 (rigid again)

**Shift map (template-based, distance to frame 0):**

* Slice 1: distance(0 vs 0) = 0 Å
* Slice 2: distance(10 vs 0) ≈ 1–2 Å (depends on frame 10)
* Slice 3: distance(20 vs 0) = 5 Å

So in the heatmaps:

* **RMSX** highlights **local flexibility in each window**.
* **Shift maps** highlight **drift from the starting structure over time**.

When you add a new metric, decide first: is it **window-based** (like RMSX) or **template-based** (like shift)?

---

## 2. What a “metric module” should look like

To keep new methods consistent with RMSX and shift maps, use this pattern:

1. **Low-level: “per-frame metric” function**

   Inputs:

   * `MDAnalysis.Universe` at the *current* frame
   * selection string (e.g. `"protein and name CA"`)

   Outputs:

   * 1D array of per-residue values
   * residue identifiers (resid + segid)

   Example skeleton:

   ```python
   def metric_for_frame(universe, selection="protein and name CA"):
       atoms    = universe.select_atoms(selection)
       residues = atoms.residues
       n_res    = len(residues)

       # Real implementation goes here:
       #   values[i] = some_function_of_frame(universe, residue_i)
       values = np.zeros(n_res, dtype=float)

       return values, residues
   ```

2. **Mid-level: “time-series driver”**

   * Loops over slices / frames.
   * Calls `metric_for_frame(...)` per slice.
   * Builds a CSV with **one row per residue** and **one column per slice**:

     ```text
     ResidueID,ChainID,slice_1,slice_2,...,slice_N
     10,A,0.12,0.15,...
     11,A,0.30,0.28,...
     ```

   Example (window or template logic swapped in as needed):

   ```python
   def compute_metric_time_series(universe,
                                  selection="protein and name CA",
                                  slice_indices=None,
                                  out_dir="metric_output"):
       os.makedirs(out_dir, exist_ok=True)

       if slice_indices is None:
           n_frames      = universe.trajectory.n_frames
           slice_indices = list(range(0, n_frames, max(1, n_frames // 10)))

       slice_names  = [f"slice_{i+1}" for i in range(len(slice_indices))]
       rows         = {}
       residues_ref = None

       for k, frame_idx in enumerate(slice_indices):
           universe.trajectory[frame_idx]
           values, residues = metric_for_frame(universe, selection=selection)

           if residues_ref is None:
               residues_ref = residues
               for res in residues_ref:
                   key = (res.resid, res.segid)
                   rows[key] = {"ResidueID": res.resid, "ChainID": res.segid}

           col_name = slice_names[k]
           for res, v in zip(residues, values):
               key = (res.resid, res.segid)
               rows[key][col_name] = float(v)

       csv_path   = os.path.join(out_dir, "metric_timeseries.csv")
       fieldnames = ["ResidueID", "ChainID"] + slice_names

       with open(csv_path, "w", newline="") as f:
           writer = csv.DictWriter(f, fieldnames=fieldnames)
           writer.writeheader()
           for key in sorted(rows.keys()):
               writer.writerow(rows[key])

       return csv_path
   ```

3. **High-level: `run_<metric>()`**

   Does the same jobs as `run_rmsx()` or `run_shift_map()`:

   * load `Universe`
   * select chain
   * choose slices / frame indices
   * call `compute_metric_time_series(...)`
   * send CSV to:

     * `file_namer(...)` for a nice filename
     * `update_all_pdb_bfactors(...)` for B-factors
     * `create_r_plot(...)` for heatmaps

   Minimal sketch:

   ```python
   def run_metric(topology_file,
                  trajectory_file,
                  out_dir=None,
                  num_slices=10,
                  selection="protein and name CA",
                  log_transform=False,
                  custom_fill_label="MyMetric"):
       u = mda.Universe(topology_file, trajectory_file)

       # decide on slice indices (window-based or template-based)
       n_frames      = len(u.trajectory)
       slice_indices = list(range(0, n_frames, max(1, n_frames // num_slices)))

       csv_raw = compute_metric_time_series(
           u,
           selection=selection,
           slice_indices=slice_indices,
           out_dir=out_dir or "metric_output",
       )

       df = pd.read_csv(csv_raw)
       if log_transform:
           slice_cols = [c for c in df.columns if c not in ("ResidueID", "ChainID")]
           df[slice_cols] = df[slice_cols].apply(np.log1p)
           df.to_csv(csv_raw, index=False)

       update_all_pdb_bfactors(csv_raw, silent=False, verbose=True)

       rmsd_csv = calculate_rmsd(u, out_dir, analysis_type="protein")
       rmsf_csv = calculate_rmsf(u, out_dir, analysis_type="protein")

       create_r_plot(
           csv_raw,
           rmsd_csv,
           rmsf_csv,
           interpolate=False,
           triple=True,
           palette="viridis",
           custom_fill_label=custom_fill_label,
       )

       return csv_raw
   ```

4. **Flipbook hook: `run_<metric>_flipbook()`**

   Wraps `all_chain_<metric>()` and then calls `run_flipbook(...)`, mirroring `run_rmsx_flipbook()` and `run_shift_flipbook()`.

---

## 3. Worked Example: Per-Residue lDDT (stored as 1 − lDDT)

lDDT is naturally **template-based**: for each frame, you compare local distances to a **reference structure** and score how many are preserved within fixed thresholds (usually 0.5, 1, 2, 4 Å). The canonical definition is “higher = better local agreement” (1 = perfect, 0 = bad).

To keep RMSX, shift maps, and lDDT visually consistent, we:

* **Compute lDDT as usual in [0, 1]**, then
* **Store `1 − lDDT`** in the CSV and PDB B-factors, so that:

  * **0** → perfect agreement / very stable
  * **1** → highly distorted / unstable
* Label the plots and legends as **“1 − lDDT (CA-only)”** so the semantics are explicit.

### 3.1. lDDT helpers (math only)

The helpers still compute **plain lDDT**, i.e. high = good. The inversion happens later.

```python
# lddt_utils.py
import numpy as np
from MDAnalysis.lib.distances import distance_array

def compute_reference_distances(universe,
                                selection="protein and name CA",
                                inclusion_radius=15.0):
    """
    Precompute CA-CA reference distances and neighbor lists.
    Assumes one CA per residue.
    """
    ca_atoms = universe.select_atoms(selection)
    residues = ca_atoms.residues
    n_res    = len(residues)

    if n_res == 0:
        raise ValueError(f"No residues found for selection: {selection}")

    coords   = ca_atoms.positions               # (N, 3)
    ref_dmat = distance_array(coords, coords)   # (N, N)

    neighbor_indices = []
    for i in range(n_res):
        mask      = (ref_dmat[i] < inclusion_radius)
        mask[i]   = False
        neighbors = np.where(mask)[0]
        neighbor_indices.append(neighbors)

    return ref_dmat, neighbor_indices, residues


def compute_lddt_per_frame(universe,
                           ref_dmat,
                           neighbor_indices,
                           selection="protein and name CA",
                           thresholds=(0.5, 1.0, 2.0, 4.0)):
    """
    Compute CA-only lDDT for the *current* frame.
    Returns scores (N,) and residue list.

    NOTE: This returns **lDDT** in [0, 1] (higher = more local agreement).
          The inversion to (1 - lDDT) happens later in the pipeline.
    """
    ca_atoms = universe.select_atoms(selection)
    residues = ca_atoms.residues
    n_res    = len(residues)

    coords = ca_atoms.positions
    dmat   = distance_array(coords, coords)  # current frame

    thr = np.asarray(thresholds, dtype=float)
    scores = np.zeros(n_res, dtype=float)

    for i in range(n_res):
        neigh = neighbor_indices[i]
        if neigh.size == 0:
            scores[i] = np.nan
            continue

        d_ref = ref_dmat[i, neigh]
        d_cur = dmat[i, neigh]
        diff  = np.abs(d_cur - d_ref)

        preserved      = diff[None, :] < thr[:, None]
        frac_preserved = preserved.sum(axis=1) / float(neigh.size)
        scores[i]      = frac_preserved.mean()

    return scores, residues
```

This module:

* Knows **nothing** about CSVs, R, or Flipbook.
* Returns **plain lDDT**; we explicitly invert it later to keep semantics consistent.

---

### 3.2. Template-based 1 − lDDT time-series (shift-style)

Now plug lDDT into the **template-based** slice driver, and invert it **once**, right after computing the scores:

```python
def compute_lddt_time_series(u,
                             ref_u,
                             selection="protein and name CA",
                             num_slices=10,
                             start_frame=0,
                             out_dir="lddt_output",
                             inclusion_radius=15.0,
                             thresholds=(0.5, 1.0, 2.0, 4.0)):
    """
    Template-based lDDT:
      - reference distances from ref_u
      - one frame per slice from u
      - output: residue × slice CSV of **1 - lDDT**,
        so that higher values = more local instability.

    The inversion (1 - lDDT) is applied once here so the rest of the
    RMSX/Flipbook pipeline can treat this like “another instability metric”.
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) Precompute reference geometry
    ref_dmat, neighbor_indices, ref_residues = compute_reference_distances(
        ref_u,
        selection=selection,
        inclusion_radius=inclusion_radius,
    )

    # 2) Decide slice frames (like shift map)
    n_frames      = len(u.trajectory) - start_frame
    slice_size    = max(1, n_frames // num_slices)
    slice_indices = [start_frame + i * slice_size for i in range(num_slices)]
    slice_names   = [f"slice_{i+1}" for i in range(num_slices)]

    rows = {}
    for k, frame_idx in enumerate(slice_indices):
        u.trajectory[frame_idx]

        # lDDT in [0, 1] (higher = better local agreement)
        lddt_scores, residues = compute_lddt_per_frame(
            u,
            ref_dmat,
            neighbor_indices,
            selection=selection,
            thresholds=thresholds,
        )

        # Invert once: instability = 1 - lDDT
        values = 1.0 - np.asarray(lddt_scores, dtype=float)

        if k == 0:
            # init rows keyed by (resid, segid)
            for res in residues:
                key = (res.resid, res.segid)
                rows[key] = {"ResidueID": res.resid, "ChainID": res.segid}

        col_name = slice_names[k]
        for res, v in zip(residues, values):
            key = (res.resid, res.segid)
            rows[key][col_name] = float(v)

    csv_path   = os.path.join(out_dir, "lddt_timeseries.csv")
    fieldnames = ["ResidueID", "ChainID"] + slice_names
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(rows.keys()):
            writer.writerow(rows[key])

    return csv_path
```

Key points:

* **Only this function** knows about the inversion.
* Everything downstream (R plots, Flipbook) just sees “higher = more unstable”, same as RMSX/shift.

---

### 3.3. High-level `run_lddt_map()` (exporting 1 − lDDT)

Finally, integrate into the RMSX pipeline with a high-level wrapper that:

* Uses `compute_lddt_time_series(...)` (which already outputs 1 − lDDT).
* Writes B-factors and calls the R heatmap script.
* Sets a clear legend label reflecting the inversion.

```python
def run_lddt_map(topology_file,
                 trajectory_file,
                 output_dir=None,
                 num_slices=10,
                 chain_sele=None,
                 log_transform=False):
    """
    High-level lDDT driver, patterned after run_shift_map().

    IMPORTANT:
      - The CSV and B-factors store **1 - lDDT (CA-only)**.
      - This keeps the visual semantics:
            higher value  -> more motion / less local agreement
        consistent with RMSX and shift maps.
    """
    u     = mda.Universe(topology_file, trajectory_file)
    ref_u = mda.Universe(topology_file)

    # simple chain selection (you can copy the full logic from run_rmsx)
    if chain_sele:
        selection = f"protein and name CA and segid {chain_sele}"
    else:
        selection = "protein and name CA"

    out_dir = output_dir or f"{Path(topology_file).stem}_lddtmap"

    # Compute **1 - lDDT** time-series
    csv_raw = compute_lddt_time_series(
        u,
        ref_u,
        selection=selection,
        num_slices=num_slices,
        start_frame=0,
        out_dir=out_dir,
    )

    df = pd.read_csv(csv_raw)
    if log_transform:
        slice_cols = [c for c in df.columns if c not in ("ResidueID", "ChainID")]
        df[slice_cols] = df[slice_cols].apply(np.log1p)
        df.to_csv(csv_raw, index=False)

    # reuse RMSX plumbing
    rmsd_csv = calculate_rmsd(u, out_dir, analysis_type="protein")
    rmsf_csv = calculate_rmsf(u, out_dir, analysis_type="protein")

    # Write **1 - lDDT** into PDB B-factors
    update_all_pdb_bfactors(csv_raw, silent=False, verbose=True)

    # Heatmap label explicitly says 1 - lDDT
    create_r_plot(
        csv_raw,
        rmsd_csv,
        rmsf_csv,
        interpolate=False,
        triple=True,
        palette="viridis",
        custom_fill_label="1 − lDDT (CA-only)",
    )

    return csv_raw
```

You then add the usual:

* `all_chain_lddt_map(...)`
* `run_lddt_flipbook(...)`

following the **exact** pattern of `all_chain_rmsx` and `run_rmsx_flipbook`, but with:

* CSVs already containing **1 − lDDT**,
* B-factors already set to **1 − lDDT**,
* and legends labelled `"1 − lDDT (CA-only)"`.

### Masking with `lDDT`

`lDDT` can use the same residue-level masking flow as RMSX and shift maps. Pass raw
MDAnalysis selection strings to `mask=` to exclude disordered regions from the
unmasked min/max calculation while still keeping them visible in the CSV, heatmap,
and Flipbook output.

```python
from rmsx.addons.lddt import run_lddt_map, run_lddt_flipbook

run_lddt_map(
    topology_file="prot.pdb",
    trajectory_file="traj.dcd",
    output_dir="lddt_masked",
    num_slices=9,
    chain_sele="A",
    mask="segid A and resid 103:150",
)

run_lddt_flipbook(
    topology_file="prot.pdb",
    trajectory_file="traj.dcd",
    output_dir="lddt_masked_flipbook",
    num_slices=9,
    mask=[
        "segid A and resid 103:150",
        "segid B and resid 12:25",
    ],
)
```

Masked `lDDT` runs behave like masked RMSX/shift runs:

* masked residues are excluded from summaries and synced color-scale estimation,
* masked rows are clipped against the unmasked floor/ceiling after any `log_transform`,
* `masked_residues.csv` is written beside each chain CSV and in the combined Flipbook directory,
* heatmaps use the same hatched overlay, and Flipbook viewers use the same transparency cues.

---

## 4. Adding Other Metrics Later

Once lDDT is in place, adding another per-residue metric is:

1. Decide: **window-based** (RMSX-style) or **template-based** (shift-style).
2. Write a **pure metric helper** module:

   * `metric_for_frame(...)` (and, if needed, `compute_reference_*()`).
3. Write a **time-series driver**:

   * `compute_<metric>_time_series(...)` → residue × slice CSV.
   * If the natural semantics are “high = good” but you want “high = unstable”,
     do the inversion in this function (like `values = 1 - score` or similar).
4. Write a **high-level driver**:

   * `run_<metric>_map(...)` + optional `run_<metric>_flipbook(...)`.
   * Reuse `file_namer`, `update_all_pdb_bfactors`, `calculate_rmsd`, `calculate_rmsf`, `create_r_plot`, and `run_flipbook`.

Because RMSX (window-based), shift map (template-based), and 1 − lDDT (template-based “instability”) all follow the same skeleton, new metrics can plug into RMSX/Flipbook with very little extra plumbing and consistent visual semantics.
