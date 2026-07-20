# rmsx/addons/lddt.py
"""
Addon: per-residue lDDT maps integrated with RMSX / Flipbook.

IMPORTANT NOTE:
    Internally we compute the standard lDDT score in [0, 1], but the values
    written to CSV/B-factors are *1 - lDDT*, so that "higher = more unstable"
    matches the RMSX convention.

This module follows the same patterns as RMSX and shift maps:

  * low-level metric helpers (compute_reference_distances, compute_lddt_per_frame)
  * a per-slice driver (process_trajectory_lddt_by_size)
  * a single-chain driver (run_lddt_map)
  * a multi-chain wrapper (all_chain_lddt_map)
  * a Flipbook entry point (run_lddt_flipbook)

Typical usage (single chain):

    from rmsx.addons.lddt import run_lddt_map

    run_lddt_map(
        topology_file="prot.pdb",
        trajectory_file="traj.dcd",
        num_slices=16,
        chain_sele="A",
    )

Flipbook-style usage (all chains, combined, tiled):

    from rmsx.addons.lddt import run_lddt_flipbook

    run_lddt_flipbook(
        topology_file="prot.pdb",
        trajectory_file="traj.dcd",
        num_slices=9,
        viewer="chimerax",
    )
"""

import os
from pathlib import Path
from typing import Iterable, Tuple, Optional, List

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from ..core import (
    apply_mask_clipping,
    build_mask_metadata,
    get_selection_string,
    has_masked_residues,
    initialize_environment,
    setup_directory,
    prepare_managed_output_dir,
    file_namer,
    finalize_masked_csv,
    calculate_rmsd,
    calculate_rmsf,
    update_all_pdb_bfactors,
    create_r_plot,
    summarize_rmsx,
    compute_global_rmsx_min_max,
    combine_pdb_files,
    write_combined_mask_metadata,
    write_mask_metadata,
)
from ..flipbook import run_flipbook


# ---------------------------------------------------------------------------
# 1. lDDT helpers: pure metric logic (no I/O)
# ---------------------------------------------------------------------------

def compute_reference_distances(
    universe: mda.Universe,
    selection: str = "protein and name CA",
    inclusion_radius: float = 15.0,
) -> Tuple[np.ndarray, List[np.ndarray], mda.core.groups.ResidueGroup]:
    """
    Precompute reference CA-CA distances + neighbor lists per residue.

    Assumes a "one CA per residue" selection in the CURRENT frame of `universe`.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe positioned at the desired reference frame.
    selection : str
        MDAnalysis selection string (usually CA-only).
    inclusion_radius : float
        Distance cutoff (Å) for neighbors contributing to lDDT.

    Returns
    -------
    ref_dists : (N, N) np.ndarray
        Pairwise CA-CA distances in the reference frame.
    neighbor_indices : list of np.ndarray
        neighbor_indices[i] is an array of residue indices j that are within
        inclusion_radius of residue i (excluding self).
    residues : MDAnalysis.core.groups.ResidueGroup
        Residues corresponding to the CA atoms.
    """
    ca_atoms = universe.select_atoms(selection)
    residues = ca_atoms.residues
    n_res = len(residues)

    if n_res == 0:
        raise ValueError(f"No residues found with selection: '{selection}'")

    ca_positions = ca_atoms.positions  # (N, 3)
    if ca_positions.shape[0] != n_res:
        raise ValueError(
            "Number of selected atoms != number of residues; "
            "this helper assumes one CA per residue."
        )

    ref_dists = distance_array(ca_positions, ca_positions)  # (N, N)

    neighbor_indices: List[np.ndarray] = []
    for i in range(n_res):
        mask = (ref_dists[i] < inclusion_radius)
        mask[i] = False  # exclude self
        neighbors = np.where(mask)[0]
        neighbor_indices.append(neighbors)

    return ref_dists, neighbor_indices, residues


def compute_lddt_per_frame(
    universe: mda.Universe,
    ref_dists: np.ndarray,
    neighbor_indices: List[np.ndarray],
    selection: str = "protein and name CA",
    thresholds: Iterable[float] = (0.5, 1.0, 2.0, 4.0),
) -> Tuple[np.ndarray, mda.core.groups.ResidueGroup]:
    """
    Compute per-residue CA-only lDDT for the CURRENT frame in `universe`.

    NOTE:
        This returns the *standard* lDDT in [0, 1]. The inversion to
        (1 - lDDT) is applied later when constructing the slice matrix.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe positioned at the frame to be scored.
    ref_dists : (N, N) np.ndarray
        Reference CA-CA distances from `compute_reference_distances`.
    neighbor_indices : list of np.ndarray
        neighbor_indices[i] are indices of neighbors for residue i.
    selection : str
        MDAnalysis selection string (must match the reference selection).
    thresholds : iterable of float
        Distance tolerances (Å) used to evaluate "preserved" contacts.

    Returns
    -------
    scores : 1D np.ndarray, length N
        lDDT-like score in [0, 1] per residue.
    residues : MDAnalysis.core.groups.ResidueGroup
        Residues corresponding to the CA atoms (same order as scores).
    """
    ca_atoms = universe.select_atoms(selection)
    residues = ca_atoms.residues
    n_res = len(residues)

    if ca_atoms.n_atoms != n_res:
        raise ValueError(
            "Number of selected atoms != number of residues; "
            "this helper assumes one CA per residue."
        )

    ca_positions = ca_atoms.positions
    dists = distance_array(ca_positions, ca_positions)

    thr = np.array(list(thresholds), dtype=float)
    scores = np.zeros(n_res, dtype=float)

    for i in range(n_res):
        neigh = neighbor_indices[i]
        if neigh.size == 0:
            scores[i] = np.nan
            continue

        d_ref = ref_dists[i, neigh]
        d_mod = dists[i, neigh]
        diff = np.abs(d_mod - d_ref)

        preserved = diff[None, :] < thr[:, None]  # (n_thr, n_neigh)
        frac_preserved = preserved.sum(axis=1) / float(neigh.size)  # per-threshold
        scores[i] = frac_preserved.mean()

    return scores, residues


# ---------------------------------------------------------------------------
# 2. Per-slice driver: lDDT over time (RMSX-style matrix, storing 1 - lDDT)
# ---------------------------------------------------------------------------

def process_trajectory_lddt_by_size(
    u: mda.Universe,
    output_dir: str,
    total_size: int,
    slice_size: int,
    chain_sele: Optional[str] = None,
    start_frame: int = 0,
    analysis_type: str = "protein",
    verbose: bool = True,
    inclusion_radius: float = 15.0,
    thresholds: Iterable[float] = (0.5, 1.0, 2.0, 4.0),
) -> Tuple[pd.DataFrame, int]:
    """
    Slice the trajectory into fixed-size chunks and compute per-residue lDDT
    for the first frame of each slice, using the first frame of the region
    as the template.

    The values stored in the returned DataFrame are *1 - lDDT* so that
    higher values mean greater local instability, matching RMSX semantics.

    Output:
      DataFrame with:
        - ResidueID
        - ChainID
        - slice_1.dcd, slice_2.dcd, ...

    Parameters
    ----------
    u : MDAnalysis.Universe
        Loaded MD trajectory (topology + coordinates).
    output_dir : str
        Directory where per-slice PDBs will be written.
    total_size : int
        Total number of frames in the analyzed range (e.g. end - start + 1).
    slice_size : int
        Number of frames per slice.
    chain_sele : str or None
        Chain ID (segid) to restrict analysis to a specific chain.
    start_frame : int
        Index of the first frame in the region of interest.
    analysis_type : str
        "protein" or "dna", passed through to get_selection_string.
    verbose : bool
        Print progress if True.
    inclusion_radius : float
        Cutoff (Å) for neighbors in lDDT calculation.
    thresholds : iterable of float
        Distance thresholds (Å) for lDDT.

    Returns
    -------
    all_data : pandas.DataFrame
        Per-residue (1 - lDDT) values per slice.
    adjusted_total_size : int
        Total number of frames actually used after truncation.
    """
    # 1) Set reference frame and selection (CA-only for metric)
    u.trajectory[start_frame]
    selection_str = get_selection_string(
        analysis_type=analysis_type,
        chain_sele=chain_sele,
        full_backbone=False,   # CA-only for lDDT
    )

    ref_dists, neighbor_indices, ref_residues = compute_reference_distances(
        u,
        selection=selection_str,
        inclusion_radius=inclusion_radius,
    )

    # 2) Normalize total_size to an integer number of slices
    adjusted_total_size = (total_size // slice_size) * slice_size
    n_slices = adjusted_total_size // slice_size

    if verbose:
        print(
            f"Computing lDDT (storing 1 - lDDT) for {n_slices} slices "
            f"from frames {start_frame} to {start_frame + adjusted_total_size - 1}."
        )

    all_data = pd.DataFrame()

    for i in range(n_slices):
        frame_idx = start_frame + i * slice_size
        u.trajectory[frame_idx]

        # lDDT scores for this frame (0 = unstable, 1 = stable)
        lddt_scores, residues = compute_lddt_per_frame(
            u,
            ref_dists,
            neighbor_indices,
            selection=selection_str,
            thresholds=thresholds,
        )

        # Invert to instability: 1 - lDDT
        instability_scores = 1.0 - lddt_scores

        col_name = f"slice_{i + 1}.dcd"
        df_slice = pd.DataFrame(
            {col_name: instability_scores},
            index=[res.resid for res in residues],
        )

        all_data = df_slice if all_data.empty else pd.concat([all_data, df_slice], axis=1)

        # Write a PDB for this slice using full backbone for Flipbook
        atoms_full = u.select_atoms(
            get_selection_string(
                analysis_type=analysis_type,
                chain_sele=chain_sele,
                full_backbone=True,   # backbone for visualization
            )
        )
        coord_path = os.path.join(output_dir, f"slice_{i + 1}_first_frame.pdb")
        with mda.Writer(coord_path, atoms_full.n_atoms, multiframe=False) as w:
            w.write(atoms_full)

        if verbose:
            print(
                f"Slice {i + 1}: lDDT computed at frame {frame_idx}, "
                f"stored as 1 - lDDT, PDB written to {coord_path}"
            )

    # 3) Add residue / chain identifiers using the reference ordering
    all_data.insert(0, "ChainID", [res.atoms[0].segid for res in ref_residues])
    all_data.insert(0, "ResidueID", [res.resid for res in ref_residues])

    return all_data, adjusted_total_size


# ---------------------------------------------------------------------------
# 3. High-level driver: run_lddt_map (single-chain, RMSX-style)
# ---------------------------------------------------------------------------

def run_lddt_map(
    topology_file: str,
    trajectory_file: str,
    output_dir: Optional[str] = None,
    num_slices: Optional[int] = None,
    slice_size: Optional[int] = None,
    rscript_executable: str = "Rscript",
    verbose: bool = True,
    interpolate: bool = True,
    triple: bool = False,
    chain_sele: Optional[str] = None,
    overwrite: bool = False,
    palette: str = "viridis",
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    make_plot: bool = True,
    analysis_type: str = "protein",
    summary_n: Optional[int] = 3,
    manual_length_ns: Optional[float] = None,
    log_transform: bool = False,
    inclusion_radius: float = 15.0,
    thresholds: Iterable[float] = (0.5, 1.0, 2.0, 4.0),
    custom_fill_label: str = "Instability\n(1 − lDDT)",
    mask=None,
    _defer_mask_clipping: bool = False,
    _allow_fully_masked: bool = False,
):
    """
    High-level driver for per-residue lDDT maps.

    NOTE:
        Values stored in the CSV/B-factors are 1 - lDDT, so that
        higher values correspond to greater local instability.

    Mirrors run_shift_map / run_rmsx:
      - chooses frames
      - calls process_trajectory_lddt_by_size(...)
      - saves CSV with file_namer(...)
      - writes B-factors
      - calls create_r_plot(...)

    Returns
    -------
    summary_tuple : (top_n_df, bottom_n_df) or None
        Top/bottom (1 − lDDT) residues, if summary_n is not None.
    """
    initialize_environment(verbose=verbose)

    # 1) Decide output directory
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_lddtmap")

    # 2) Chain selection (using topology only)
    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)

    if chain_sele is None:
        chain_info = {}
        for chain in chain_ids:
            chain_atoms = u_top.select_atoms(f"segid {chain}")
            chain_info[chain] = len(chain_atoms.residues)
        if verbose:
            print("Available chains and their lengths (in residues):")
            for chain, length in chain_info.items():
                print(f"  Chain {chain}: {length} residues")

        chain_list = ", ".join(
            [f"{chain} ({length} residues)" for chain, length in chain_info.items()]
        )
        selected_chain = input(
            "Please enter the chain ID you would like to analyze from the following options:\n"
            f"{chain_list}\nChain ID: "
        ).strip()
        if selected_chain not in chain_ids:
            raise RuntimeError("Selected chain is not available in the topology.")
        chain_sele = selected_chain
    else:
        if chain_sele not in chain_ids:
            raise RuntimeError(f"Chain '{chain_sele}' is not available in the topology.")

    # 3) Prepare chain-specific output directory
    chain_output_dir = os.path.join(output_dir, f"chain_{chain_sele}_lddtmap")
    setup_directory(chain_output_dir, overwrite=overwrite, verbose=verbose)

    # 4) Full trajectory universe
    u = mda.Universe(topology_file, trajectory_file)

    if end_frame is None:
        end_frame = len(u.trajectory) - 1
    if end_frame < start_frame:
        raise ValueError("end_frame must be >= start_frame.")
    if end_frame >= len(u.trajectory):
        end_frame = len(u.trajectory) - 1

    used_frames_count = end_frame - start_frame + 1
    if used_frames_count < 1:
        raise ValueError("No frames available in the specified range.")

    if verbose:
        print("Starting lDDT map analysis (storing 1 − lDDT)...")

    # 5) Determine slice_size from num_slices if needed
    if slice_size is None:
        if num_slices is not None:
            slice_size = used_frames_count // num_slices
            if verbose:
                print(f"Calculated slice_size = {slice_size} (from num_slices = {num_slices})")
        else:
            raise RuntimeError("Either slice_size or num_slices must be specified.")

    # 6) Compute lDDT matrix per slice (stored as 1 - lDDT)
    all_data, adjusted_total_size = process_trajectory_lddt_by_size(
        u,
        chain_output_dir,
        used_frames_count,
        slice_size,
        chain_sele=chain_sele,
        start_frame=start_frame,
        analysis_type=analysis_type,
        verbose=verbose,
        inclusion_radius=inclusion_radius,
        thresholds=thresholds,
    )

    # 7) Optional log-transform (like RMSX)
    if log_transform:
        slice_cols = [c for c in all_data.columns if c not in ("ResidueID", "ChainID")]
        all_data[slice_cols] = all_data[slice_cols].apply(np.log1p)
        new_names = {col: col.replace(".dcd", "_log.dcd") for col in slice_cols}
        all_data.rename(columns=new_names, inplace=True)

    # 8) Save CSV using the same naming logic (but prefix 'lddt')
    lddt_csv = file_namer(
        chain_output_dir,
        trajectory_file,
        "csv",
        u=u,
        frames_used=adjusted_total_size,
        manual_length_ns=manual_length_ns,
        prefix="lddt",
    )

    mask_metadata = build_mask_metadata(
        universe=u,
        selection_str=get_selection_string(
            analysis_type=analysis_type,
            chain_sele=chain_sele,
            full_backbone=False,
        ),
        data_frame=all_data,
        mask=mask,
        verbose=verbose,
    )

    if _defer_mask_clipping:
        clipped_data = all_data
    else:
        clipped_data, _ = apply_mask_clipping(
            all_data,
            mask_metadata,
            allow_fully_masked=_allow_fully_masked,
        )

    clipped_data.to_csv(lddt_csv, index=False)
    write_mask_metadata(mask_metadata, lddt_csv)
    if verbose:
        print(f"lDDT (stored as 1 − lDDT) data saved to {lddt_csv}")

    # 9) RMSD / RMSF for context
    rmsd_csv = calculate_rmsd(
        u,
        chain_output_dir,
        chain_sele=chain_sele,
        start_frame=start_frame,
        end_frame=(start_frame + adjusted_total_size - 1),
        analysis_type=analysis_type,
        verbose=verbose,
    )
    rmsf_csv = calculate_rmsf(
        u,
        chain_output_dir,
        chain_sele=chain_sele,
        start_frame=start_frame,
        end_frame=(start_frame + adjusted_total_size - 1),
        analysis_type=analysis_type,
        verbose=verbose,
    )

    # 10) Write (1 − lDDT) values into PDB B-factors
    update_all_pdb_bfactors(lddt_csv, silent=(not verbose), verbose=verbose)

    # 11) Heatmap via existing R pipeline
    if make_plot:
        if verbose:
            print("Generating (1 − lDDT) plots...")
        create_r_plot(
            lddt_csv,
            rmsd_csv,
            rmsf_csv,
            rscript_executable=rscript_executable,
            interpolate=interpolate,
            triple=triple,
            palette=palette,
            verbose=verbose,
            log_transform=log_transform,
            custom_fill_label=custom_fill_label,
        )
    else:
        if verbose:
            print("Skipping plot generation in run_lddt_map() because make_plot=False.")

    # 12) Optional top/bottom summary (reuses summarize_rmsx)
    summary_tuple = None
    if summary_n is not None and isinstance(summary_n, int):
        if verbose:
            print(
                f"\nNow summarizing the top {summary_n} and "
                f"bottom {summary_n} (1 − lDDT) values..."
            )
        top_n_df, bottom_n_df = summarize_rmsx(lddt_csv, n=summary_n, print_output=verbose)
        summary_tuple = (top_n_df, bottom_n_df)

    return summary_tuple


# ---------------------------------------------------------------------------
# 4. Multi-chain wrapper: all_chain_lddt_map (pattern of all_chain_rmsx)
# ---------------------------------------------------------------------------

def all_chain_lddt_map(
    topology_file: str,
    trajectory_file: str,
    output_dir: Optional[str] = None,
    num_slices: Optional[int] = None,
    slice_size: Optional[int] = None,
    rscript_executable: str = "Rscript",
    verbose: bool = True,
    interpolate: bool = True,
    triple: bool = False,
    overwrite: bool = False,
    palette: str = "viridis",
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    sync_color_scale: bool = True,
    analysis_type: str = "protein",
    manual_length_ns: Optional[float] = None,
    summary_n: Optional[int] = 3,
    log_transform: bool = False,
    inclusion_radius: float = 15.0,
    thresholds: Iterable[float] = (0.5, 1.0, 2.0, 4.0),
    custom_fill_label: str = "Instability\n(1 − lDDT)",
    mask=None,
) -> str:
    """
    Perform lDDT map analysis for all chains in the topology file.

    Values stored are 1 − lDDT, so that higher values indicate more local
    instability, matching RMSX's "high = floppy" convention.

    Follows the pattern of all_chain_rmsx / all_chain_shift_map:

      * runs run_lddt_map(...) for each segid
      * optionally syncs the color scale across chains
      * combines chain-specific PDBs into a single 'combined' directory

    Returns
    -------
    combined_dir : str
        Directory containing combined PDBs (for Flipbook).
    """
    mask_is_active = has_masked_residues(mask)
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_lddtmap")

    prepare_managed_output_dir(
        output_dir,
        overwrite=overwrite,
        verbose=verbose,
        topology_file=topology_file,
        trajectory_file=trajectory_file,
    )

    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)

    combined_output_dirs = []
    csv_paths: List[str] = []

    for chain in chain_ids:
        if verbose:
            print(f"\nAnalyzing Chain {chain}...")

        chain_make_plot = (not sync_color_scale) and (not mask_is_active)
        _summary_tuple = run_lddt_map(
            topology_file=topology_file,
            trajectory_file=trajectory_file,
            output_dir=output_dir,
            num_slices=num_slices,
            slice_size=slice_size,
            rscript_executable=rscript_executable,
            verbose=verbose,
            interpolate=interpolate,
            triple=triple,
            chain_sele=chain,
            overwrite=False,
            palette=palette,
            start_frame=start_frame,
            end_frame=end_frame,
            make_plot=chain_make_plot,
            analysis_type=analysis_type,
            summary_n=None if mask_is_active else summary_n,
            manual_length_ns=manual_length_ns,
            log_transform=log_transform,
            inclusion_radius=inclusion_radius,
            thresholds=thresholds,
            custom_fill_label=custom_fill_label,
            mask=mask,
            _defer_mask_clipping=mask_is_active,
            _allow_fully_masked=mask_is_active,
        )

        chain_output_dir = os.path.join(output_dir, f"chain_{chain}_lddtmap")
        combined_output_dirs.append(chain_output_dir)

        possible_csv = list(Path(chain_output_dir).glob("lddt_*.csv"))
        if possible_csv:
            csv_paths.append(str(possible_csv[0]))

    if mask_is_active:
        if verbose:
            print("\nComputing global (1 − lDDT) min and max across all unmasked residues...")
        global_min, global_max = compute_global_rmsx_min_max(csv_paths)
        if verbose:
            print(f"Global (1 − lDDT) range = [{global_min:.3f}, {global_max:.3f}]")

        for csv_path in csv_paths:
            finalize_masked_csv(csv_path, (global_min, global_max), verbose=verbose)

        if summary_n is not None and isinstance(summary_n, int):
            for csv_path in csv_paths:
                chain_label = Path(csv_path).parent.name
                if verbose:
                    print(f"\nFinal masked summary for {chain_label}...")
                summarize_rmsx(csv_path, n=summary_n, print_output=verbose)

    # Combine PDBs across chains (for multi-chain Flipbook)
    if len(chain_ids) > 1:
        combined_dir = os.path.join(output_dir, "combined")
        if verbose:
            print("\nCombining PDB files from all chains...")
        combine_pdb_files(combined_output_dirs, combined_dir, verbose=verbose)
        if mask_is_active:
            write_combined_mask_metadata(combined_output_dirs, combined_dir)
        if verbose:
            print("Combined lDDT analysis completed for all chains.")
    else:
        single_chain_id = chain_ids[0]
        combined_dir = os.path.join(output_dir, f"chain_{single_chain_id}_lddtmap")
        if mask_is_active:
            write_combined_mask_metadata(combined_output_dirs, combined_dir)
        if verbose:
            print(f"Single-chain lDDT analysis completed. Using directory: {combined_dir}")

    # Optional global color scale sync
    if mask_is_active and csv_paths:
        if verbose:
            if sync_color_scale:
                print("Generating final masked (1 − lDDT) plots with a fixed color scale...")
            else:
                print("Generating final masked (1 − lDDT) plots...")

        for csv_path in csv_paths:
            csv_dir = Path(csv_path).parent
            rmsd_csv = csv_dir / "rmsd.csv"
            rmsf_csv = csv_dir / "rmsf.csv"
            create_r_plot(
                rmsx_csv=str(csv_path),
                rmsd_csv=str(rmsd_csv),
                rmsf_csv=str(rmsf_csv),
                rscript_executable=rscript_executable,
                interpolate=interpolate,
                triple=triple,
                palette=palette,
                min_val=global_min if sync_color_scale else None,
                max_val=global_max if sync_color_scale else None,
                verbose=verbose,
                log_transform=log_transform,
                custom_fill_label=custom_fill_label,
            )

    elif sync_color_scale and csv_paths:
        if verbose:
            print("\nComputing global (1 − lDDT) min and max across all chains...")
        global_min, global_max = compute_global_rmsx_min_max(csv_paths)
        if verbose:
            print(f"Global (1 − lDDT) range = [{global_min:.3f}, {global_max:.3f}]")
            print("Generating final (1 − lDDT) plots with a fixed color scale...")

        for csv_path in csv_paths:
            csv_dir = Path(csv_path).parent
            rmsd_csv = csv_dir / "rmsd.csv"
            rmsf_csv = csv_dir / "rmsf.csv"
            create_r_plot(
                rmsx_csv=str(csv_path),
                rmsd_csv=str(rmsd_csv),
                rmsf_csv=str(rmsf_csv),
                rscript_executable=rscript_executable,
                interpolate=interpolate,
                triple=triple,
                palette=palette,
                min_val=global_min,
                max_val=global_max,
                verbose=verbose,
                log_transform=log_transform,
                custom_fill_label=custom_fill_label,
            )

    return combined_dir


# ---------------------------------------------------------------------------
# 5. Flipbook hook: run_lddt_flipbook (pattern of run_rmsx_flipbook)
# ---------------------------------------------------------------------------

def run_lddt_flipbook(
    topology_file: str,
    trajectory_file: str,
    output_dir: Optional[str] = None,
    num_slices: int = 9,
    slice_size: Optional[int] = None,
    rscript_executable: str = "Rscript",
    verbose: bool = True,
    interpolate: bool = True,
    triple: bool = False,
    overwrite: bool = False,
    palette: str = "viridis",
    spacingFactor: Optional[str] = None,
    start_frame: int = 0,
    end_frame: Optional[int] = None,
    analysis_type: str = "protein",
    manual_length_ns: Optional[float] = None,
    summary_n: Optional[int] = 3,
    sync_color_scale: bool = True,
    flipbook_min_bfactor: Optional[float] = None,
    flipbook_max_bfactor: Optional[float] = None,
    log_transform: bool = False,
    inclusion_radius: float = 15.0,
    thresholds: Iterable[float] = (0.5, 1.0, 2.0, 4.0),
    custom_fill_label: str = "Instability\n(1 − lDDT)",
    viewer: str = "chimerax",
    extra_commands: Optional[Iterable[str]] = None,
    mask=None,
    molstar_output=None,
    molstar_manifest_output=None,
    molstar_asset_mode: str = "cdn",
    molstar_height: int = 720,
    molstar_camera_mode: str = "orthographic",
    molstar_orientation: str = "portrait",
):
    """
    Run lDDT analysis and generate a Flipbook visualization, syncing the
    color scale across all chains by default (mirrors run_rmsx_flipbook).

    Values passed to Flipbook via B-factors are 1 − lDDT, so thicker,
    hotter worms correspond to locally less accurate regions.

    Parameters
    ----------
    viewer : str
        "molstar", "chimerax", or "vmd", passed through to run_flipbook.
    extra_commands : iterable of str or None
        Extra commands to append in the viewer script (e.g., ChimeraX commands).

    Returns
    -------
    combined_dir : str
        Directory used as input to Flipbook (combined PDBs).
    """
    combined_dir = all_chain_lddt_map(
        topology_file=topology_file,
        trajectory_file=trajectory_file,
        output_dir=output_dir,
        num_slices=num_slices,
        slice_size=slice_size,
        rscript_executable=rscript_executable,
        verbose=verbose,
        interpolate=interpolate,
        triple=triple,
        overwrite=overwrite,
        palette=palette,
        start_frame=start_frame,
        end_frame=end_frame,
        sync_color_scale=sync_color_scale,
        analysis_type=analysis_type,
        manual_length_ns=manual_length_ns,
        summary_n=summary_n,
        log_transform=log_transform,
        inclusion_radius=inclusion_radius,
        thresholds=thresholds,
        custom_fill_label=custom_fill_label,
        mask=mask,
    )

    result = run_flipbook(
        directory=combined_dir,
        palette=palette,
        min_bfactor=flipbook_min_bfactor,
        max_bfactor=flipbook_max_bfactor,
        spacingFactor=spacingFactor,
        extra_commands=extra_commands,
        viewer=viewer,
        molstar_output=molstar_output,
        molstar_manifest_output=molstar_manifest_output,
        molstar_asset_mode=molstar_asset_mode,
        molstar_height=molstar_height,
        molstar_camera_mode=molstar_camera_mode,
        molstar_orientation=molstar_orientation,
    )

    if verbose:
        print("Full (1 − lDDT) Flipbook analysis completed successfully.")

    if str(viewer).lower() == "molstar":
        return result
    return combined_dir
