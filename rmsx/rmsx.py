import warnings  # Import warnings first to suppress specific warnings before other imports

# ---------------------------------------------------------------------
#                  Original notes/documentation
# ---------------------------------------------------------------------
# The code below provides a framework for analyzing MD trajectories using MDAnalysis.
# It slices trajectories, computes RMSF/RMSD, updates PDB B-factors, and can optionally
# create RMSX plots in R and produce a 3D 'FlipBook'-style visualization.
#
# We have:
#  - Maintained all original comments, notes, and code structure.
#  - Added `summarize_rmsx(...)` to compute and optionally print top & bottom residues.
#  - Added `summary_n` (default=3) to `run_rmsx(...)` so it returns top/bottom RMSX values if desired.
#  - Added `manual_length_ns` parameter in `file_namer` so we can override the auto-calculated simulation length.
#  - Passed a `verbose` parameter around to control console output (only prints if `verbose=True`).
#  - Ensured you can override the default selection string (e.g. for DNA) by passing `atom_selection`.
#
# If you want to see the top 5 & bottom 5 RMSX residues, call run_rmsx(..., summary_n=5).
# If you want a manual simulation length in ns, specify manual_length_ns=100 (for example).
# If you don’t want verbose output, set verbose=False.


# ---------------------------------------------------------------------
#                  Notes / Documentation
# ---------------------------------------------------------------------
# This script slices trajectories, computes RMSF/RMSD, updates PDB B-factors,
# can create RMSX plots via an R script, produce a 3D "FlipBook"-style visualization,
# and summarize top/bottom RMSX residues.
#
# Key features:
# 1. Analysis type parameter (analysis_type="protein" or "nucleic") for selecting
#    either alpha-carbons or nucleic acid atoms.
# 2. summary_n: if set to an integer, automatically summarize the top/bottom RMSX values.
# 3. manual_length_ns: if set, we override the auto-calculated simulation length
#    (based on dt * frames) and use the user-specified length in the output filename.
# 4. verbose: only prints progress messages if True (default). If False, minimal output.
# ---------------------------------------------------------------------

warnings.filterwarnings(
    "ignore",  # Completely ignore the warning
    message="Found no information for attr: 'formalcharges' Using default value of '0'",
    category=UserWarning,
    module="MDAnalysis.coordinates.PDB"
)

import os
import re
import sys
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.rms import RMSF
import subprocess
import glob
import shutil
import numpy as np
import plotly.graph_objects as go

from IPython.display import Image, display
from contextlib import redirect_stdout
from pathlib import Path

# If you have the flipbook module in the same directory structure, update import as needed
from .flipbook import run_flipbook


def get_selection_string(analysis_type="protein", chain_sele=None):
    """
    Return the MDAnalysis selection string for a given analysis type.
    Defaults to protein (CA atoms).
    For DNA, you can choose e.g. "nucleic and name P" or "nucleic and name C4'".

    analysis_type: str
        Either "protein" (default) or "dna" (or another defined string).
    chain_sele: str
        If provided, we append 'and segid <chain_sele>'.
    """
    if analysis_type.lower() == "dna":
        # Example: select phosphate atoms in nucleic acid
        base_selection = "nucleic and name P"
    else:
        # Default to protein alpha carbons
        base_selection = "protein and name CA"

    if chain_sele:
        base_selection += f" and segid {chain_sele}"

    return base_selection


def summarize_rmsx(rmsx_csv, n=3, print_output=True):
    """
    Read the RMSX CSV file and return a summary of the top N and bottom N
    RMSX values across all slices (columns).

    Parameters
    ----------
    rmsx_csv : str
        Path to the RMSX CSV file.
    n : int, optional
        Number of top/bottom entries to show (default=3).
    print_output : bool, optional
        If True (default), prints the summary to the console.

    Returns
    -------
    top_n_df : pd.DataFrame
        DataFrame of the top N rows, with columns: [ResidueID, ChainID, TimeSlice, RMSX].
    bottom_n_df : pd.DataFrame
        DataFrame of the bottom N rows, with columns: [ResidueID, ChainID, TimeSlice, RMSX].
    """
    df = pd.read_csv(rmsx_csv)

    # Identify numeric columns except 'ResidueID' and 'ChainID'
    skip_cols = ['ResidueID', 'ChainID']
    numeric_cols = [c for c in df.columns if c not in skip_cols]

    # Melt the dataframe so each slice (column) becomes a row
    melted = df.melt(
        id_vars=skip_cols,
        value_vars=numeric_cols,
        var_name='TimeSlice',
        value_name='RMSX'
    )

    # Sort to find top/bottom n
    top_n_df = melted.nlargest(n, 'RMSX')
    bottom_n_df = melted.nsmallest(n, 'RMSX')

    if print_output:
        print(f"\n===== Top {n} RMSX values across all slices =====")
        for _, row in top_n_df.iterrows():
            print(
                f"TimeSlice: {row['TimeSlice']}, "
                f"ChainID: {row['ChainID']}, "
                f"ResidueID: {row['ResidueID']}, "
                f"RMSX: {row['RMSX']:.3f}"
            )

        print(f"\n===== Bottom {n} RMSX values across all slices =====")
        for _, row in bottom_n_df.iterrows():
            print(
                f"TimeSlice: {row['TimeSlice']}, "
                f"ChainID: {row['ChainID']}, "
                f"ResidueID: {row['ResidueID']}, "
                f"RMSX: {row['RMSX']:.3f}"
            )
        print()

    return top_n_df, bottom_n_df


def initialize_environment(verbose=True):
    """Print the Python executable path and current working directory, if verbose."""
    if verbose:
        print("Python Executable:", sys.executable)
        print("Current Working Directory:", os.getcwd())


def check_directory(output_dir, overwrite=False, verbose=True):
    """
    Check if the output directory exists and handle overwriting based on the overwrite flag.
    Prints if 'verbose' is True; otherwise remains silent (except for essential user input).
    """
    if os.path.exists(output_dir):
        if overwrite:
            if verbose:
                print(f"Overwriting existing directory: {output_dir}")
            return True
        else:
            # User input is essential, so we do not wrap it in 'verbose'
            response = input(
                f"The directory '{output_dir}' already exists. Overwrite? (y/n): "
            )
            return response.strip().lower() == 'y'
    return True  # Proceed if the directory does not exist


def setup_directory(output_dir, overwrite=False, verbose=True):
    """
    Set up or clear the output directory based on user preference.
    """
    if check_directory(output_dir, overwrite=overwrite, verbose=verbose):
        if os.path.exists(output_dir):
            # Clear existing directory content
            for file in os.listdir(output_dir):
                file_path = os.path.join(output_dir, file)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    if verbose:
                        print(f'Failed to delete {file_path}. Reason: {e}')
        else:
            os.makedirs(output_dir)
        if verbose:
            print(f"The directory '{output_dir}' is ready for use.")
    else:
        if verbose:
            print(
                f"Process terminated by user. The directory '{output_dir}' will not be overwritten."
            )
        raise RuntimeError("User chose not to overwrite the existing directory.")


def process_trajectory_slices_by_size(
    u,
    output_dir,
    total_size,
    slice_size,
    chain_sele=None,
    start_frame=0,
    analysis_type="protein",
    verbose=True
):
    """
    Slice the trajectory based on a fixed number of frames per slice, ensuring all slices have the same size.
    Truncate excess frames if total_size is not divisible by slice_size.
    """
    adjusted_total_size = (total_size // slice_size) * slice_size
    excess_frames = total_size - adjusted_total_size

    if excess_frames > 0:
        if verbose:
            print(f"Truncating {excess_frames} excess frame(s). Original size: {total_size}, Updated: {adjusted_total_size}")
        total_size = adjusted_total_size
    else:
        if verbose:
            print(f"No truncation needed. Total size: {total_size} frames")

    n_slices = total_size // slice_size
    if n_slices == 0:
        raise ValueError("Slice size is larger than the total number of frames available in the chosen range.")

    if verbose:
        print(f"Processing frames {start_frame} to {start_frame + total_size - 1} of the trajectory.")
        print(f"Number of slices: {n_slices}")

    selection_str = get_selection_string(analysis_type=analysis_type, chain_sele=chain_sele)
    atoms_to_analyze = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    # for i in range(n_slices):
    #     slice_start = start_frame + i * slice_size
    #     slice_end = slice_start + slice_size  # stop is exclusive in RMSF.run
    #     # Move to the start frame and write the first-frame PDB
    #     u.trajectory[slice_start]
    #     coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
    #     with mda.Writer(coord_path, atoms_to_analyze.n_atoms, multiframe=False) as coord_writer:
    #         coord_writer.write(atoms_to_analyze)
    #     if verbose:
    #         print(f"First frame of slice {i + 1} written to {coord_path}")
    #
    #     rmsf_calc = RMSF(atoms_to_analyze)
    #     rmsf_calc.run(start=slice_start, stop=slice_end)

    #seeing if this fixes the issue with psf loading of files.
    for i in range(n_slices):
        slice_start = start_frame + i * slice_size
        slice_end = slice_start + slice_size  # stop is exclusive in RMSF.run

        # Update to the desired frame
        u.trajectory[slice_start]

        # Re-run the selection to get the current frame’s coordinates
        atoms_to_analyze = u.select_atoms(selection_str)
        if len(atoms_to_analyze) == 0:
            raise ValueError(
                f"Selection returned empty at frame {slice_start}. Check your selection or coordinate file.")

        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, atoms_to_analyze.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(atoms_to_analyze)
        if verbose:
            print(f"First frame of slice {i + 1} written to {coord_path}")

        rmsf_calc = RMSF(atoms_to_analyze)
        rmsf_calc.run(start=slice_start, stop=slice_end)

        df = pd.DataFrame(
            {f"slice_{i + 1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in atoms_to_analyze.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        if verbose:
            print(f"Slice {i + 1}: RMSF computed for frames {slice_start} to {slice_end - 1} ({slice_size} frames)")

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(0, 'ChainID', [res.atoms[0].segid for res in atoms_to_analyze.residues])
        all_data.insert(0, 'ResidueID', [res.resid for res in atoms_to_analyze.residues])

    return all_data, adjusted_total_size


def process_trajectory_slices_by_num(
    u,
    output_dir,
    total_size,
    num_slices,
    chain_sele=None,
    start_frame=0,
    analysis_type="protein",
    verbose=True
):
    """
    Slice the trajectory into a specific number of slices, ensuring all slices have the same size.
    Truncate excess frames if total_size is not divisible by num_slices.
    """
    adjusted_total_size = (total_size // num_slices) * num_slices
    excess_frames = total_size - adjusted_total_size

    if excess_frames > 0:
        if verbose:
            print(f"Truncating {excess_frames} excess frame(s). Original size: {total_size}, Updated: {adjusted_total_size}")
        total_size = adjusted_total_size
    else:
        if verbose:
            print(f"No truncation needed. Total size: {total_size} frames")

    if verbose:
        print(f"Processing frames {start_frame} to {start_frame + total_size - 1} of the trajectory.")
        print(f"Number of slices: {num_slices}")

    base_size = total_size // num_slices
    slice_sizes = [base_size] * num_slices  # All slices have the same size

    selection_str = get_selection_string(analysis_type=analysis_type, chain_sele=chain_sele)
    atoms_to_analyze = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    current_start = start_frame
    for i, size in enumerate(slice_sizes):
        slice_start = current_start
        slice_end = slice_start + size
        current_start += size

        u.trajectory[slice_start]
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, atoms_to_analyze.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(atoms_to_analyze)
        if verbose:
            print(f"First frame of slice {i + 1} written to {coord_path}")

        rmsf_calc = RMSF(atoms_to_analyze)
        rmsf_calc.run(start=slice_start, stop=slice_end)

        df = pd.DataFrame(
            {f"slice_{i + 1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in atoms_to_analyze.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        if verbose:
            print(f"Slice {i + 1}: RMSF computed for frames {slice_start} to {slice_end - 1} ({size} frames)")

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(0, 'ChainID', [res.atoms[0].segid for res in atoms_to_analyze.residues])
        all_data.insert(0, 'ResidueID', [res.resid for res in atoms_to_analyze.residues])

    return all_data, adjusted_total_size


def analyze_trajectory(output_dir, chain_sele=None):
    """(Not used in this approach.)"""
    return pd.DataFrame()


def file_namer(
    output_dir,
    example_file,
    out_file_type="csv",
    prefix="rmsx",
    u=None,
    frames_used=None,
    manual_length_ns=None
):
    """
    Generate a filename based on simulation parameters.

    If manual_length_ns is provided, override the auto-calculated length.
    """
    if manual_length_ns is not None:
        # If user manually specifies length, just use it
        simulation_length_ns = manual_length_ns
    else:
        simulation_length_fs = extract_simulation_length(u, frames_used=frames_used)
        simulation_length_ns = simulation_length_fs / 1e6  # Convert fs to ns

    sim_name = os.path.basename(example_file)
    sim_name = os.path.splitext(sim_name)[0]

    if simulation_length_ns == 0:
        decimals = 3
    elif simulation_length_ns < 0.001:
        decimals = 6
    elif simulation_length_ns < 0.01:
        decimals = 5
    else:
        decimals = 3

    output_filename = f'{prefix}_{sim_name}_{simulation_length_ns:.{decimals}f}_ns.{out_file_type}'
    output_filepath = os.path.join(output_dir, output_filename)
    return output_filepath


def extract_simulation_length(u, frames_used=None):
    """
    Extract simulation length in fs based on frames_used.
    If frames_used is None, use entire trajectory.
    """
    timestep_fs = u.trajectory.dt * 1000.0  # ps -> fs
    if frames_used is None:
        frames_used = len(u.trajectory)
    total_time_fs = timestep_fs * frames_used
    return total_time_fs


def calculate_rmsd(
    u,
    output_dir,
    selection='protein and name CA',
    chain_sele=None,
    start_frame=0,
    end_frame=None,
    analysis_type="protein",
    verbose=True
):
    """
    Calculate RMSD for the subset of the trajectory defined by start_frame and end_frame.
    """
    selection_str = get_selection_string(analysis_type=analysis_type, chain_sele=chain_sele)
    atoms_to_analyze = u.select_atoms(selection_str)

    rmsd_analysis = rms.RMSD(atoms_to_analyze)
    stop_frame = None
    if end_frame is not None:
        stop_frame = end_frame + 1
    rmsd_analysis.run(start=start_frame, stop=stop_frame)

    columns = ['Frame', 'Time', 'RMSD']
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd, columns=columns)

    if output_dir:
        rmsd_output_filepath = os.path.join(output_dir, 'rmsd.csv')
        rmsd_df.to_csv(rmsd_output_filepath, index=False)
        if verbose:
            print(f"RMSD data saved to {rmsd_output_filepath}")
        return rmsd_output_filepath


def calculate_rmsf(
    u,
    output_dir=None,
    selection='protein and name CA',
    chain_sele=None,
    start_frame=0,
    end_frame=None,
    analysis_type="protein",
    verbose=True
):
    """
    Calculate RMSF for the subset of the trajectory defined by start_frame and end_frame.
    """
    u.trajectory[start_frame]

    selection_str = get_selection_string(analysis_type=analysis_type, chain_sele=chain_sele)
    atoms_to_analyze = u.select_atoms(selection_str)

    rmsf_analysis = RMSF(atoms_to_analyze)
    stop_frame = None
    if end_frame is not None:
        stop_frame = end_frame + 1
    rmsf_analysis.run(start=start_frame, stop=stop_frame)
    rmsf_values = rmsf_analysis.results.rmsf

    num_residues = len(atoms_to_analyze.residues)
    if verbose:
        print(f"Number of residues: {num_residues}")
    rmsf_list = rmsf_values.tolist()
    if len(rmsf_list) != num_residues:
        raise ValueError(
            f"Length of rmsf_list ({len(rmsf_list)}) != number of residues ({num_residues})"
        )

    rmsf_whole_traj = pd.DataFrame(
        {
            'ResidueID': [residue.resid for residue in atoms_to_analyze.residues],
            'RMSF': rmsf_list,
        }
    )

    if output_dir:
        rmsf_output_filepath = os.path.join(output_dir, 'rmsf.csv')
        rmsf_whole_traj.to_csv(rmsf_output_filepath, index=False)
        if verbose:
            print(f"RMSF data saved to {rmsf_output_filepath}")
        return rmsf_output_filepath


def create_r_plot(
    rmsx_csv,
    rmsd_csv,
    rmsf_csv,
    rscript_executable='Rscript',
    interpolate=False,
    triple=False,
    palette="plasma",
    min_val=None,
    max_val=None,
    verbose=True
):
    """
    Run the R script to generate RMSX plots and display the first image.
    """
    interpolate_str = 'TRUE' if interpolate else 'FALSE'
    triple_str = 'TRUE' if triple else 'FALSE'

    try:
        try:
            current_dir = Path(__file__).parent.resolve()
        except NameError:
            current_dir = Path.cwd().resolve()

        r_script_path = current_dir / 'r_scripts' / 'plot_rmsx.R'
        if not r_script_path.is_file():
            if verbose:
                print(f"Error: R script not found at {r_script_path}.")
            return

        if verbose:
            print(f"Found R script at {r_script_path}.")

        cmd = [
            rscript_executable,
            str(r_script_path),
            rmsx_csv,
            rmsd_csv if rmsd_csv else "",
            rmsf_csv if rmsf_csv else "",
            interpolate_str,
            triple_str,
            palette,
            str(min_val) if min_val is not None else "",
            str(max_val) if max_val is not None else ""
        ]

        if verbose:
            print("Running R script command:")
            print(" ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            if verbose:
                print("R script executed successfully.")
        else:
            if verbose:
                print("R script execution failed.")
                print(f"Error Output:\n{result.stderr}")
            return

    except FileNotFoundError:
        if verbose:
            print(
                f"Error: Rscript executable not found: {rscript_executable}. "
                "Please ensure R is installed and 'Rscript' is in your PATH."
            )
        return
    except Exception as e:
        if verbose:
            print(f"An unexpected error occurred: {e}")
        return

    # Display the first PNG found
    try:
        rmsx_path = Path(rmsx_csv).resolve()
        directory = rmsx_path.parent

        if not directory.is_dir():
            if verbose:
                print(f"Error: The directory does not exist: {directory}")
            return

        image_files = list(directory.glob('*.png'))
        if image_files:
            first_image = image_files[0]
            if verbose:
                print(f"Displaying image: {first_image}")
            display(Image(filename=str(first_image)))
        else:
            if verbose:
                print("No PNG files found in the specified directory.")

    except Exception as e:
        if verbose:
            print(f"An error occurred while searching for PNG files: {e}")


def update_pdb_bfactor(coord_file, rmsx_df, silent=False, verbose=True):
    """
    Update the B-factor field in a PDB file with RMSX values.
    If 'silent' is True, do not print. If 'silent' is False, print only if 'verbose' is True.
    """
    coord_file_base_name = os.path.basename(coord_file)
    slice_number = int(coord_file_base_name.split('_')[1])
    rmsx_column = f'slice_{slice_number}.dcd'

    if rmsx_column not in rmsx_df.columns:
        if not silent and verbose:
            print(f"Error: {rmsx_column} not found in RMSX DataFrame.")
        return

    with open(coord_file, 'r') as pdb:
        pdb_lines = pdb.readlines()

    updated_pdb_lines = []
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            residue_number = int(line[22:26].strip())
            if residue_number in rmsx_df['ResidueID'].values:
                rmsx_value = rmsx_df.loc[
                    rmsx_df['ResidueID'] == residue_number, rmsx_column
                ].values[0]
                new_bfactor = f"{rmsx_value:6.2f}"
                updated_line = line[:60] + new_bfactor + line[66:]
                updated_pdb_lines.append(updated_line)
            else:
                updated_pdb_lines.append(line)
        else:
            updated_pdb_lines.append(line)

    with open(coord_file, 'w') as pdb:
        pdb.writelines(updated_pdb_lines)
    if not silent and verbose:
        print(f"Updated PDB B-factors in {coord_file}")


def load_coord_files(folder_path):
    """
    Load and sort PDB coordinate files from a folder.
    """
    def extract_number(filename):
        match = re.search(r'\d+', filename)
        return int(match.group()) if match else float('inf')

    coord_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    sorted_coord_files = sorted(coord_files, key=extract_number)
    return [os.path.join(folder_path, f) for f in sorted_coord_files]


def update_all_pdb_bfactors(rmsx_csv, silent, verbose=True):
    """
    Update B-factors for all PDB files in the directory based on RMSX data.
    """
    rmsx_df = pd.read_csv(rmsx_csv)
    coord_folder = os.path.dirname(rmsx_csv)
    coord_files = load_coord_files(coord_folder)
    for coord_file in coord_files:
        update_pdb_bfactor(coord_file, rmsx_df, silent=silent, verbose=verbose)


def combine_pdb_files(chain_dirs, combined_dir, silent=False, verbose=True):
    """
    Combine PDB files from multiple chains into a single PDB file per slice.
    """
    os.makedirs(combined_dir, exist_ok=True)
    pdb_files = [f for f in os.listdir(chain_dirs[0]) if f.endswith('.pdb')]
    for pdb_file in pdb_files:
        combined_content = []
        missing_files = False
        for chain_dir in chain_dirs:
            pdb_file_path = os.path.join(chain_dir, pdb_file)
            if os.path.exists(pdb_file_path):
                with open(pdb_file_path, 'r') as file:
                    combined_content.extend(file.readlines())
            else:
                if verbose:
                    print(f"File {pdb_file} not found in {chain_dir}. Skipping.")
                missing_files = True
                break
        if not missing_files:
            combined_pdb_file = os.path.join(combined_dir, pdb_file)
            with open(combined_pdb_file, 'w') as combined_file:
                combined_file.writelines(combined_content)
            if not silent and verbose:
                print(f"Combined {pdb_file} from all chains into {combined_pdb_file}")


def find_and_combine_pdb_files(output_dir, verbose=True):
    """
    Find all chain-specific RMSX directories and combine their PDB files.
    """
    directories_with_dir = [
        os.path.join(output_dir, d)
        for d in os.listdir(output_dir)
        if d.endswith('_rmsx')
           and not d.startswith("combined")
           and os.path.isdir(os.path.join(output_dir, d))
    ]
    if verbose:
        print("Directories to combine:", directories_with_dir)
    combined_dir = os.path.join(output_dir, "combined")
    combine_pdb_files(directories_with_dir, combined_dir, verbose=verbose)
    if verbose:
        print(f"Combined PDB files are saved in '{combined_dir}'.")


def compute_global_rmsx_min_max(csv_paths):
    """
    Read all provided RMSX CSV files and compute the overall min and max RMSX values.
    """
    import numpy as np

    all_vals = []
    for path in csv_paths:
        df = pd.read_csv(path)
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col not in ('ResidueID', 'ChainID'):
                all_vals.append(df[col].values)

    if not all_vals:
        return 0.0, 1.0

    all_vals = np.concatenate(all_vals)
    global_min = float(all_vals.min())
    global_max = float(all_vals.max())
    return global_min, global_max


def run_rmsx(
    topology_file,
    trajectory_file,
    output_dir=None,
    num_slices=None,
    slice_size=None,
    rscript_executable='Rscript',
    verbose=True,
    interpolate=True,
    triple=False,
    chain_sele=None,
    overwrite=False,
    palette="viridis",
    start_frame=0,
    end_frame=None,
    make_plot=True,
    analysis_type="protein",
    summary_n=3,
    manual_length_ns=None
):
    """
    Run the RMSX analysis on a specified trajectory range.

    Parameters
    ----------
    - topology_file : str
    - trajectory_file : str
    - output_dir : str or None
    - num_slices : int or None
    - slice_size : int or None
    - rscript_executable : str
    - verbose : bool
    - interpolate : bool
    - triple : bool
    - chain_sele : str or None
    - overwrite : bool
    - palette : str
    - start_frame : int
    - end_frame : int or None
    - make_plot : bool
    - analysis_type : str ("protein" or "dna")
    - summary_n : int or None
        If provided, we automatically compute and return top/bottom RMSX data.
    - manual_length_ns : float or None
        If provided, overrides the auto-calculated simulation length for naming output files.

    Returns
    -------
    summary_tuple : (pd.DataFrame, pd.DataFrame) or None
        If summary_n is not None, returns (top_n_df, bottom_n_df). Otherwise None.
    """
    initialize_environment(verbose=verbose)

    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")

    # Load a minimal universe for chain info
    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)
    chain_info = {}
    for chain in chain_ids:
        chain_atoms = u_top.select_atoms(f'segid {chain}')
        num_residues = len(chain_atoms.residues)
        chain_info[chain] = num_residues

    if chain_sele is None:
        if verbose:
            print("Available chains and their lengths (in residues):")
            for chain, length in chain_info.items():
                print(f"Chain {chain}: {length} residues")
            chain_list = ", ".join(
                [f"{chain} ({length} residues)" for chain, length in chain_info.items()]
            )
        selected_chain = input(
            f"Please enter the chain ID you would like to analyze from the following options:\n{chain_list}\nChain ID: "
        ).strip()
        if selected_chain not in chain_ids:
            if verbose:
                print(f"Chain '{selected_chain}' is not available in the topology file.")
            raise RuntimeError("Selected chain is not available in the topology.")
        chain_sele = selected_chain
    else:
        if chain_sele not in chain_ids:
            if verbose:
                print(f"Chain '{chain_sele}' is not available in the topology file.")
            raise RuntimeError("Selected chain is not available in the topology.")

    base_name = os.path.splitext(os.path.basename(topology_file))[0]
    chain_output_dir = os.path.join(output_dir, f"chain_{chain_sele}_rmsx")
    setup_directory(chain_output_dir, overwrite=overwrite, verbose=verbose)

    # Create the main universe with trajectory
    u = mda.Universe(topology_file, trajectory_file)

    # If end_frame not specified, use last frame
    if end_frame is None:
        end_frame = len(u.trajectory) - 1

    if end_frame < start_frame:
        if verbose:
            print("Error: end_frame must be >= start_frame.")
        raise ValueError("end_frame must be >= start_frame.")
    if end_frame >= len(u.trajectory):
        end_frame = len(u.trajectory) - 1

    used_frames_count = end_frame - start_frame + 1
    if used_frames_count < 1:
        if verbose:
            print("No frames available in the specified range.")
        raise ValueError("No frames available in the specified range.")

    if verbose:
        print("Starting analysis...")

    # Slicing approach
    if num_slices is not None:
        if verbose:
            print(f"Using the slicing method with num_slices={num_slices}")
        all_data, adjusted_total_size = process_trajectory_slices_by_num(
            u,
            chain_output_dir,
            used_frames_count,
            num_slices,
            chain_sele=chain_sele,
            start_frame=start_frame,
            analysis_type=analysis_type,
            verbose=verbose
        )
    elif slice_size is not None:
        if verbose:
            print(f"Using the slicing method with slice_size={slice_size}")
        all_data, adjusted_total_size = process_trajectory_slices_by_size(
            u,
            chain_output_dir,
            used_frames_count,
            slice_size,
            chain_sele=chain_sele,
            start_frame=start_frame,
            analysis_type=analysis_type,
            verbose=verbose
        )
        # We haven't saved the CSV yet, so do that now:
        # The actual saving is below, after we create the filename
    else:
        if verbose:
            print("Error: You must specify either num_slices or slice_size.")
        raise RuntimeError("No slicing method specified.")

    # Create the RMSX CSV filename
    rmsx_csv = file_namer(
        chain_output_dir,
        trajectory_file,
        "csv",
        u=u,
        frames_used=adjusted_total_size,
        manual_length_ns=manual_length_ns
    )
    # Save the combined data
    all_data.to_csv(rmsx_csv, index=False)
    if verbose:
        print(f"RMSX data saved to {rmsx_csv}")

    # RMSD
    rmsd_csv = calculate_rmsd(
        u,
        chain_output_dir,
        chain_sele=chain_sele,
        start_frame=start_frame,
        end_frame=(start_frame + adjusted_total_size - 1),
        analysis_type=analysis_type,
        verbose=verbose
    )
    # RMSF
    rmsf_csv = calculate_rmsf(
        u,
        chain_output_dir,
        chain_sele=chain_sele,
        start_frame=start_frame,
        end_frame=(start_frame + adjusted_total_size - 1),
        analysis_type=analysis_type,
        verbose=verbose
    )

    # Update PDB B-factors
    update_all_pdb_bfactors(rmsx_csv, silent=(not verbose), verbose=verbose)

    # Optional plot
    if make_plot:
        if verbose:
            print("Generating plots...")
        create_r_plot(
            rmsx_csv,
            rmsd_csv,
            rmsf_csv,
            rscript_executable=rscript_executable,
            interpolate=interpolate,
            triple=triple,
            palette=palette,
            verbose=verbose
        )
    else:
        if verbose:
            print("Skipping plot generation in run_rmsx() because make_plot=False.")

    # Optional summary
    summary_tuple = None
    if summary_n is not None and isinstance(summary_n, int):
        if verbose:
            print(f"\nNow summarizing the top {summary_n} and bottom {summary_n} RMSX values...")
        top_n_df, bottom_n_df = summarize_rmsx(rmsx_csv, n=summary_n, print_output=verbose)
        summary_tuple = (top_n_df, bottom_n_df)

    return summary_tuple


def all_chain_rmsx(
    topology_file,
    trajectory_file,
    output_dir=None,
    num_slices=None,
    slice_size=None,
    rscript_executable='Rscript',
    verbose=True,
    interpolate=True,
    triple=False,
    overwrite=False,
    palette='viridis',
    start_frame=0,
    end_frame=None,
    sync_color_scale=False,
    analysis_type="protein",
    manual_length_ns=None,
    summary_n=3
):
    """
    Perform RMSX analysis for all chains in the topology file.

    If sync_color_scale=True, we skip immediate plotting in run_rmsx
    and do a global min/max pass after analyzing all chains.
    """
    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)

    combined_output_dirs = []
    csv_paths = []

    for chain in chain_ids:
        if verbose:
            print(f"\nAnalyzing Chain {chain}...")

        chain_make_plot = not sync_color_scale
        # Run RMSX for each chain
        _summary_tuple = run_rmsx(
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
            overwrite=overwrite,
            palette=palette,
            start_frame=start_frame,
            end_frame=end_frame,
            make_plot=chain_make_plot,
            analysis_type=analysis_type,
            summary_n=summary_n,
            manual_length_ns=manual_length_ns
        )

        chain_output_dir = os.path.join(output_dir, f"chain_{chain}_rmsx")
        combined_output_dirs.append(chain_output_dir)

        possible_csv = list(Path(chain_output_dir).glob("rmsx_*.csv"))
        if possible_csv:
            csv_paths.append(str(possible_csv[0]))

    if len(chain_ids) > 1:
        combined_dir = os.path.join(output_dir, "combined")
        if verbose:
            print("\nCombining PDB files from all chains...")
        combine_pdb_files(combined_output_dirs, combined_dir, verbose=verbose)
        if verbose:
            print("Combined RMSX analysis completed for all chains.")
    else:
        single_chain_id = chain_ids[0]
        combined_dir = os.path.join(output_dir, f"chain_{single_chain_id}_rmsx")
        if verbose:
            print(f"Single-chain analysis completed. Using directory: {combined_dir}")

    # If sync_color_scale, do a global pass for min/max and re-plot
    if sync_color_scale and csv_paths:
        if verbose:
            print("\nComputing global RMSX min and max across all chains...")
        global_min, global_max = compute_global_rmsx_min_max(csv_paths)
        if verbose:
            print(f"Global RMSX range = [{global_min:.3f}, {global_max:.3f}]")
            print("Generating final plots with a fixed color scale...")

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
                verbose=verbose
            )

    return combined_dir


def run_rmsx_flipbook(
    topology_file,
    trajectory_file,
    output_dir=None,
    num_slices=9,
    slice_size=None,
    rscript_executable='Rscript',
    verbose=True,
    interpolate=True,
    triple=False,
    overwrite=False,
    palette='viridis',
    spacingFactor="1",
    start_frame=0,
    end_frame=None,
    analysis_type="protein and name CA",
    manual_length_ns=None,
    summary_n=3,
    flipbook_min_bfactor=None,
    flipbook_max_bfactor=None
):
    """
    Run RMSX analysis and generate a FlipBook visualization, syncing the color scale
    across all chains by default.
    """
    combined_dir = all_chain_rmsx(
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
        sync_color_scale=True,
        analysis_type=analysis_type,
        manual_length_ns=manual_length_ns,
        summary_n=summary_n
    )

    run_flipbook(
        directory=combined_dir,
        palette=palette,
        min_bfactor=flipbook_min_bfactor,
        max_bfactor=flipbook_max_bfactor,
        spacingFactor=spacingFactor
    )

    if verbose:
        print("Full analysis including FlipBook visualization completed successfully.")


# ---------------------------------------------------------------------
# Example usage (commented out):
#
# run_rmsx_flipbook(
#     topology_file="myprotein.pdb",
#     trajectory_file="mytrajectory.dcd",
#     output_dir="my_output",
#     num_slices=10,       # or slice_size=100
#     rscript_executable="Rscript",
#     verbose=True,
#     interpolate=False,
#     triple=True,
#     overwrite=True,
#     palette="magma",
#     flipbook_min_bfactor=0,
#     flipbook_max_bfactor=5,
#     spacingFactor="0.5",
#     start_frame=0,
#     end_frame=None,
#     analysis_type="protein",
#     manual_length_ns=100,  # Force naming as 100 ns
#     summary_n=5            # Summarize top/bottom 5
# )
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Example usage for single-chain RMSX analysis with uniform slice sizes
# (Commented out; update paths/parameters to your environment)
#
# traj_file_single = "/path/to/your/single_chain_trajectory.dcd"
# pdb_file_single = "/path/to/your/single_chain_protein.pdb"
# output_dir_single = "/path/to/your/output_directory_single_chain"
#
# run_rmsx_flipbook(
#     topology_file=pdb_file_single,
#     trajectory_file=traj_file_single,
#     output_dir=output_dir_single,
#     num_slices=16,
#     slice_size=None,
#     rscript_executable='Rscript',
#     verbose=True,
#     interpolate=False,
#     triple=True,
#     overwrite=True,
#     palette="magma",
#     spacingFactor="0.5",
#     start_frame=0,
#     end_frame=None
# )
#
#
# Example usage for multi-chain RMSX analysis:
# traj_file_multi = "/path/to/your/multi_chain_trajectory.xtc"
# pdb_file_multi = "/path/to/your/multi_chain_protein.pdb"
# output_dir_multi = "/path/to/your/output_directory_multi_chain"
#
# run_rmsx_flipbook(
#     topology_file=pdb_file_multi,
#     trajectory_file=traj_file_multi,
#     output_dir=output_dir_multi,
#     num_slices=12,
#     slice_size=None,
#     rscript_executable='Rscript',
#     verbose=True,
#     interpolate=False,
#     triple=True,
#     overwrite=True,
#     palette="magma",
#     spacingFactor="1",
#     start_frame=0,
#     end_frame=None
# )
#
# Additional comments and to-do items:
#   - Test with more MD file types to ensure there's no issue reading the simulation duration
#   - Could add 'hover over' to show time slice in the flipbook, etc.
#   - Let user easily cut the simulation shorter by specifying frames
#   - Possibly create a function that auto-finds output data for 3D plotting or flipbook
#   - Consider applying log to RMSX columns for log-scale color mapping
#   - Option for flipbook to handle single chain more seamlessly
# ---------------------------------------------------------------------
