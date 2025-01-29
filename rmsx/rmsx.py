import warnings  # Import warnings first to suppress specific warnings before other imports

# Suppress the specific MDAnalysis warning about 'formalcharges'
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
from IPython.display import Image, display
from contextlib import redirect_stdout
import glob
import shutil
import numpy as np
import plotly.graph_objects as go

from pathlib import Path
from .flipbook import run_flipbook  # Use a relative import


def initialize_environment():
    """Print the Python executable path and current working directory."""
    print("Python Executable:", sys.executable)
    print("Current Working Directory:", os.getcwd())


def check_directory(output_dir, overwrite=False):
    """Check if the output directory exists and handle overwriting based on the overwrite flag."""
    if os.path.exists(output_dir):
        if overwrite:
            return True
        else:
            response = input(
                f"The directory '{output_dir}' already exists. Do you want to overwrite it? (y/n): "
            )
            return response.strip().lower() == 'y'
    return True  # Proceed if the directory does not exist


def setup_directory(output_dir, overwrite=False):
    """Set up or clear the output directory based on user preference."""
    if check_directory(output_dir, overwrite=overwrite):
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
                    print(f'Failed to delete {file_path}. Reason: {e}')
        else:
            os.makedirs(output_dir)
        print(f"The directory '{output_dir}' is ready for use.")
    else:
        print(
            f"Process terminated by user. The directory '{output_dir}' will not be overwritten."
        )
        raise RuntimeError("User chose not to overwrite the existing directory.")


def process_trajectory_slices_by_size(u, output_dir, total_size, slice_size, chain_sele=None, start_frame=0):
    """
    Slice the trajectory based on a fixed number of frames per slice, ensuring all slices have the same size.
    Truncate excess frames if total_size is not divisible by slice_size.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str)
    - total_size (int): Number of frames to process
    - slice_size (int)
    - chain_sele (str)
    - start_frame (int): The first frame to process

    Returns:
    A tuple containing:
    - all_data (pd.DataFrame): RMSF values for each slice.
    - adjusted_total_size (int): Number of frames used after truncation.
    """
    adjusted_total_size = (total_size // slice_size) * slice_size
    excess_frames = total_size - adjusted_total_size

    if excess_frames > 0:
        print(f"Truncating {excess_frames} excess frame(s) to make total_size divisible by slice_size.")
        print(f"Original simulation size: {total_size} frames")
        print(f"Updated simulation size: {adjusted_total_size} frames")
        print(f"Frames lost: {excess_frames} frames")
        total_size = adjusted_total_size
    else:
        print(f"No truncation needed. Total size: {total_size} frames")

    n_slices = total_size // slice_size
    if n_slices == 0:
        raise ValueError("Slice size is larger than the total number of frames available in the chosen range.")

    print(f"Processing frames {start_frame} to {start_frame + total_size - 1} of the trajectory.")
    print(f"Number of slices: {n_slices}")

    selection_str = "protein and name CA"
    if chain_sele:
        selection_str += f" and segid {chain_sele}"
    protein = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    for i in range(n_slices):
        slice_start = start_frame + i * slice_size
        slice_end = slice_start + slice_size  # stop is exclusive in RMSF.run
        # Move to the start frame and write the first-frame PDB
        u.trajectory[slice_start]
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, protein.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(protein)
        print(f"First frame of slice {i + 1} written to {coord_path}")

        rmsf_calc = RMSF(protein)
        rmsf_calc.run(start=slice_start, stop=slice_end)

        df = pd.DataFrame(
            {f"slice_{i + 1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in protein.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        print(f'Slice {i + 1}: RMSF computed for frames {slice_start} to {slice_end - 1} (total {slice_size} frames)')

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(0, 'ChainID', [res.atoms[0].segid for res in protein.residues])
        all_data.insert(0, 'ResidueID', [res.resid for res in protein.residues])

    return all_data, adjusted_total_size  # Return adjusted_total_size


def process_trajectory_slices_by_num(u, output_dir, total_size, num_slices, chain_sele=None, start_frame=0):
    """
    Slice the trajectory into a specific number of slices, ensuring all slices have the same size.
    Truncate excess frames if total_size is not divisible by num_slices.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str)
    - total_size (int): Number of frames to process
    - num_slices (int)
    - chain_sele (str)
    - start_frame (int): The first frame to process

    Returns:
    A tuple containing:
    - all_data (pd.DataFrame): RMSF values for each slice.
    - adjusted_total_size (int): Number of frames used after truncation.
    """
    adjusted_total_size = (total_size // num_slices) * num_slices
    excess_frames = total_size - adjusted_total_size

    if excess_frames > 0:
        print(f"Truncating {excess_frames} excess frame(s) to make total_size divisible by num_slices.")
        print(f"Original simulation size: {total_size} frames")
        print(f"Updated simulation size: {adjusted_total_size} frames")
        print(f"Frames lost: {excess_frames} frames")
        total_size = adjusted_total_size
    else:
        print(f"No truncation needed. Total size: {total_size} frames")

    print(f"Processing frames {start_frame} to {start_frame + total_size - 1} of the trajectory.")
    print(f"Number of slices: {num_slices}")

    base_size = total_size // num_slices
    slice_sizes = [base_size] * num_slices  # All slices have the same size

    selection_str = "protein and name CA"
    if chain_sele:
        selection_str += f" and segid {chain_sele}"
    protein = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    current_start = start_frame
    for i, size in enumerate(slice_sizes):
        slice_start = current_start
        slice_end = slice_start + size
        current_start += size

        u.trajectory[slice_start]
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, protein.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(protein)
        print(f"First frame of slice {i + 1} written to {coord_path}")

        rmsf_calc = RMSF(protein)
        rmsf_calc.run(start=slice_start, stop=slice_end)

        df = pd.DataFrame(
            {f"slice_{i + 1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in protein.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        print(f'Slice {i + 1}: RMSF computed for frames {slice_start} to {slice_end - 1} (total {size} frames)')

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(0, 'ChainID', [res.atoms[0].segid for res in protein.residues])
        all_data.insert(0, 'ResidueID', [res.resid for res in protein.residues])

    return all_data, adjusted_total_size  # Return adjusted_total_size


def process_rmsx_by_slice_size(u, output_dir, start_frame, end_frame, slice_size=5, chain_sele=None,
                               trajectory_file=None):
    """
    Process RMSx by specifying the number of frames per slice within a given frame range.
    Ensures all slices have the same number of frames by truncating excess frames.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str)
    - start_frame (int)
    - end_frame (int)
    - slice_size (int)
    - chain_sele (str)
    - trajectory_file (str)

    Returns:
    None
    """
    total_size = end_frame - start_frame + 1
    all_data, adjusted_total_size = process_trajectory_slices_by_size(
        u, output_dir, total_size, slice_size, chain_sele, start_frame
    )

    # Save the combined data
    save_data(all_data, output_dir, trajectory_file, u, frames_used=adjusted_total_size)


def process_rmsx_by_num_slices(u, output_dir, start_frame, end_frame, num_slices=5, chain_sele=None,
                               trajectory_file=None):
    """
    Process RMSx by specifying the number of slices over a given frame range.
    Ensures all slices have the same number of frames by truncating excess frames.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str)
    - start_frame (int)
    - end_frame (int)
    - num_slices (int)
    - chain_sele (str)
    - trajectory_file (str)

    Returns:
    A tuple containing:
    - all_data (pd.DataFrame)
    - adjusted_total_size (int)
    """
    total_size = end_frame - start_frame + 1
    all_data, adjusted_total_size = process_trajectory_slices_by_num(
        u, output_dir, total_size, num_slices, chain_sele, start_frame
    )

    # Save the combined data
    save_data(all_data, output_dir, trajectory_file, u, frames_used=adjusted_total_size)
    return all_data, adjusted_total_size  # Ensure this return statement exists


def analyze_trajectory(output_dir, chain_sele=None):
    """Not used in this approach."""
    return pd.DataFrame()


def save_data(all_data, output_dir, trajectory_file, u, frames_used):
    """
    Save the RMSX data to a CSV file with a filename reflecting the simulation length.

    Parameters:
    - all_data (pd.DataFrame): RMSX data.
    - output_dir (str): Directory to save the CSV.
    - trajectory_file (str): Name of the trajectory file.
    - u (MDAnalysis.Universe): Universe object for timestep information.
    - frames_used (int): Number of frames used in the analysis.
    """
    output_filepath = file_namer(output_dir, trajectory_file, "csv", u=u, frames_used=frames_used)
    all_data.to_csv(output_filepath, index=False)
    print(f"RMSX data saved to {output_filepath}")


def file_namer(output_dir, example_file, out_file_type="csv", prefix="rmsx", u=None, frames_used=None):
    """
    Generate a filename based on simulation parameters.

    Parameters:
    - output_dir (str)
    - example_file (str)
    - out_file_type (str)
    - prefix (str)
    - u (MDAnalysis.Universe)
    - frames_used (int)

    Returns:
    - output_filepath (str)
    """
    simulation_length_fs = extract_simulation_length(u, frames_used=frames_used)
    sim_name = os.path.basename(example_file)
    sim_name = os.path.splitext(sim_name)[0]
    simulation_length_ns = simulation_length_fs / 1e6  # Convert fs to ns
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

    Parameters:
    - u (MDAnalysis.Universe)
    - frames_used (int)

    Returns:
    - total_time_fs (float)
    """
    timestep_fs = u.trajectory.dt * 1000.0  # ps to fs
    if frames_used is None:
        frames_used = len(u.trajectory)
    total_time_fs = timestep_fs * frames_used
    return total_time_fs


def calculate_rmsd(u, output_dir, selection='protein and name CA', chain_sele=None, start_frame=0, end_frame=None):
    """
    Calculate RMSD for the subset of the trajectory defined by start_frame and end_frame.

    If end_frame is None, uses the entire trajectory from start_frame onward.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str): Directory to save the RMSD CSV.
    - selection (str): Atom selection string.
    - chain_sele (str): Specific chain selection.
    - start_frame (int)
    - end_frame (int)

    Returns:
    - rmsd_output_filepath (str): Path to the RMSD CSV file.
    """
    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)

    rmsd_analysis = rms.RMSD(protein)
    stop_frame = None
    if end_frame is not None:
        # RMSD.run uses stop as exclusive
        stop_frame = end_frame + 1
    rmsd_analysis.run(start=start_frame, stop=stop_frame)

    columns = ['Frame', 'Time', 'RMSD']
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd, columns=columns)

    if output_dir:
        rmsd_output_filepath = os.path.join(output_dir, 'rmsd.csv')
        rmsd_df.to_csv(rmsd_output_filepath, index=False)
        print(f"RMSD data saved to {rmsd_output_filepath}")
        return rmsd_output_filepath


def calculate_rmsf(u, output_dir=None, selection='protein and name CA', chain_sele=None, start_frame=0, end_frame=None):
    """
    Calculate RMSF for the subset of the trajectory defined by start_frame and end_frame.

    If end_frame is None, uses all frames from start_frame onward.

    Parameters:
    - u (MDAnalysis.Universe)
    - output_dir (str): Directory to save the RMSF CSV.
    - selection (str): Atom selection string.
    - chain_sele (str): Specific chain selection.
    - start_frame (int)
    - end_frame (int)

    Returns:
    - rmsf_output_filepath (str): Path to the RMSF CSV file.
    """
    u.trajectory[start_frame]

    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)

    rmsf_analysis = RMSF(protein)
    stop_frame = None
    if end_frame is not None:
        stop_frame = end_frame + 1  # exclusive
    rmsf_analysis.run(start=start_frame, stop=stop_frame)
    rmsf_values = rmsf_analysis.results.rmsf

    num_residues = len(protein.residues)
    print(f"Number of residues: {num_residues}")
    rmsf_list = rmsf_values.tolist()
    if len(rmsf_list) != num_residues:
        raise ValueError(
            f"Length of rmsf_list ({len(rmsf_list)}) does not match number of residues ({num_residues})"
        )

    rmsf_whole_traj = pd.DataFrame(
        {
            'ResidueID': [residue.resid for residue in protein.residues],
            'RMSF': rmsf_list,
        }
    )

    if output_dir:
        rmsf_output_filepath = os.path.join(output_dir, 'rmsf.csv')
        rmsf_whole_traj.to_csv(rmsf_output_filepath, index=False)
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
        min_val=None,  # <-- new: optional min scale
        max_val=None  # <-- new: optional max scale
):
    """
    Run the R script to generate RMSX plots and display the first image.

    Parameters:
    - rmsx_csv (str): Path to the RMSX CSV file.
    - rmsd_csv (str): Path to the RMSD CSV file.
    - rmsf_csv (str): Path to the RMSF CSV file.
    - rscript_executable (str): Path to the Rscript executable.
    - interpolate (bool): Whether to interpolate RMSX values.
    - triple (bool): Whether to generate a triple plot (RMSX, RMSD, RMSF).
    - palette (str): Color palette for the plots.
    - min_val (float or None): If provided, sets a global min for the color scale.
    - max_val (float or None): If provided, sets a global max for the color scale.
    """
    interpolate_str = 'TRUE' if interpolate else 'FALSE'
    triple_str = 'TRUE' if triple else 'FALSE'

    try:
        # Attempt to locate the R script
        try:
            current_dir = Path(__file__).parent.resolve()
        except NameError:
            current_dir = Path.cwd().resolve()

        r_script_path = current_dir / 'r_scripts' / 'plot_rmsx.R'

        if not r_script_path.is_file():
            print(f"Error: R script not found at {r_script_path}.")
            return

        print(f"Found R script at {r_script_path}.")

        # Build the command
        cmd = [
            rscript_executable,
            str(r_script_path),
            rmsx_csv,
            rmsd_csv if rmsd_csv else "",
            rmsf_csv if rmsf_csv else "",
            interpolate_str,
            triple_str,
            palette,
            str(min_val) if min_val is not None else "",  # pass empty if not set
            str(max_val) if max_val is not None else ""
        ]

        print("Running R script command:")
        print(" ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            print("R script executed successfully.")
        else:
            print("R script execution failed.")
            print(f"Error Output:\n{result.stderr}")
            return

    except FileNotFoundError:
        print(
            f"Error: Rscript executable not found: {rscript_executable}. "
            "Please ensure R is installed and 'Rscript' is in your PATH."
        )
        return
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return

    # Display the first PNG found in the output directory
    try:
        rmsx_path = Path(rmsx_csv).resolve()
        directory = rmsx_path.parent

        if not directory.is_dir():
            print(f"Error: The directory does not exist: {directory}")
            return

        image_files = list(directory.glob('*.png'))

        if image_files:
            first_image = image_files[0]
            print(f"Displaying image: {first_image}")
            display(Image(filename=str(first_image)))
        else:
            print("No PNG files found in the specified directory.")

    except Exception as e:
        print(f"An error occurred while searching for PNG files: {e}")


def update_pdb_bfactor(coord_file, rmsx_df, silent=False):
    """
    Update the B-factor field in a PDB file with RMSX values.

    Parameters:
    - coord_file (str): Path to the PDB file.
    - rmsx_df (pd.DataFrame): DataFrame containing RMSX values.
    - silent (bool): If True, suppress print statements.
    """
    coord_file_base_name = os.path.basename(coord_file)
    slice_number = int(
        coord_file_base_name.split('_')[1]
    )
    rmsx_column = f'slice_{slice_number}.dcd'

    if rmsx_column not in rmsx_df.columns:
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
    if not silent:
        print(f"Original coordinate file {coord_file} has been updated with new B-factors.")


def load_coord_files(folder_path):
    """
    Load and sort PDB coordinate files from a folder.

    Parameters:
    - folder_path (str): Path to the folder containing PDB files.

    Returns:
    - sorted_coord_files (list): Sorted list of PDB file paths.
    """

    def extract_number(filename):
        match = re.search(r'\d+', filename)
        return int(match.group()) if match else float('inf')

    coord_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    sorted_coord_files = sorted(coord_files, key=extract_number)
    return [os.path.join(folder_path, f) for f in sorted_coord_files]


def update_all_pdb_bfactors(rmsx_csv, silent):
    """
    Update B-factors for all PDB files in the directory based on RMSX data.

    Parameters:
    - rmsx_csv (str): Path to the RMSX CSV file.
    - silent (bool): If True, suppress print statements.
    """
    rmsx_df = pd.read_csv(rmsx_csv)
    coord_folder = os.path.dirname(rmsx_csv)
    coord_files = load_coord_files(coord_folder)
    for coord_file in coord_files:
        update_pdb_bfactor(coord_file, rmsx_df, silent)


def combine_pdb_files(chain_dirs, combined_dir, silent=False):
    """
    Combine PDB files from multiple chains into a single PDB file per slice.

    Parameters:
    - chain_dirs (list): List of directories containing chain-specific PDB files.
    - combined_dir (str): Directory to save combined PDB files.
    - silent (bool): If True, suppress print statements.
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
                print(f"File {pdb_file} not found in {chain_dir}. Skipping.")
                missing_files = True
                break
        if not missing_files:
            combined_pdb_file = os.path.join(combined_dir, pdb_file)
            with open(combined_pdb_file, 'w') as combined_file:
                combined_file.writelines(combined_content)
            if not silent:
                print(f"Combined {pdb_file} from all chains into {combined_pdb_file}")


def find_and_combine_pdb_files(output_dir):
    """
    Find all chain-specific RMSX directories and combine their PDB files.

    Parameters:
    - output_dir (str): Main output directory containing chain-specific directories.
    """
    directories_with_dir = [
        os.path.join(output_dir, d)
        for d in os.listdir(output_dir)
        if d.endswith('_rmsx')
           and not d.startswith("combined")
           and os.path.isdir(os.path.join(output_dir, d))
    ]
    print("Directories to combine:", directories_with_dir)
    combined_dir = os.path.join(output_dir, "combined")
    combine_pdb_files(directories_with_dir, combined_dir)
    print(f"Combined PDB files are saved in '{combined_dir}'.")


def compute_global_rmsx_min_max(csv_paths):
    """
    Read all provided RMSX CSV files and compute the overall min and max RMSX values.

    Parameters:
    - csv_paths (list of str): Paths to RMSX CSV files.

    Returns:
    - (float, float): The global minimum and maximum across all numeric RMSX columns.
    """
    import numpy as np

    all_vals = []
    for path in csv_paths:
        df = pd.read_csv(path)
        # Identify all numeric columns except 'ResidueID' and 'ChainID'
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        # Gather those values
        for col in numeric_cols:
            if col not in ('ResidueID', 'ChainID'):
                all_vals.append(df[col].values)

    if not all_vals:
        # Fallback if no numeric columns found or empty data
        return 0.0, 1.0

    all_vals = np.concatenate(all_vals)
    global_min = float(all_vals.min())
    global_max = float(all_vals.max())
    return global_min, global_max


def run_rmsx(
        topology_file,
        trajectory_file,
        output_dir=None,
        num_slices=None,  # Optional: specify either num_slices or slice_size
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
        make_plot=True  # <-- NEW: if False, skip immediate R plotting
):
    """
    Run the RMSX analysis on a specified trajectory range.

    Parameters:
    - topology_file (str): Path to the topology file.
    - trajectory_file (str): Path to the trajectory file.
    - output_dir (str): Output directory for results.
    - num_slices (int): Number of slices to divide frames into.
    - slice_size (int): Size of each slice (in frames).
    - rscript_executable (str): Path to Rscript executable.
    - verbose (bool): Enable detailed logging.
    - interpolate (bool): Enable interpolation in plots.
    - triple (bool): Enable triple plot (RMSX, RMSD, RMSF).
    - chain_sele (str): Specific chain selection. If None, prompts user.
    - overwrite (bool): Overwrite existing output directory.
    - palette (str): Color palette for plots.
    - start_frame (int): The frame at which to start the analysis.
    - end_frame (int): The last frame to include (inclusive).
    - make_plot (bool): If True (default), run the R script to produce plots immediately.
                       If False, skip plotting so the data can be combined later.
    """
    initialize_environment()
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
            print(f"Chain '{selected_chain}' is not available in the topology file.")
            raise RuntimeError("Selected chain is not available in the topology.")
        chain_sele = selected_chain
    else:
        if chain_sele not in chain_ids:
            print(f"Chain '{chain_sele}' is not available in the topology file.")
            raise RuntimeError("Selected chain is not available in the topology.")

    base_name = os.path.splitext(os.path.basename(topology_file))[0]
    chain_output_dir = os.path.join(output_dir, f"chain_{chain_sele}_rmsx")
    setup_directory(chain_output_dir, overwrite=overwrite)

    # Create the main universe once with trajectory
    u = mda.Universe(topology_file, trajectory_file)

    # If end_frame not specified, use last frame
    if end_frame is None:
        end_frame = len(u.trajectory) - 1

    if end_frame < start_frame:
        raise ValueError("end_frame must be greater than or equal to start_frame.")
    if end_frame >= len(u.trajectory):
        end_frame = len(u.trajectory) - 1  # Cap at the last available frame

    used_frames_count = end_frame - start_frame + 1

    if used_frames_count < 1:
        raise ValueError("No frames available in the specified range.")

    if verbose:
        print("Starting analysis...")

    if num_slices is not None:
        if verbose:
            print(f"Using the slicing method with num_slices={num_slices}")
        all_data, adjusted_total_size = process_rmsx_by_num_slices(
            u, chain_output_dir, start_frame, end_frame, num_slices=num_slices, chain_sele=chain_sele,
            trajectory_file=trajectory_file
        )
    elif slice_size is not None:
        if verbose:
            print(f"Using the slicing method with slice_size={slice_size}")
        # Here, process_rmsx_by_slice_size does not return all_data, adjusted_total_size,
        # but we can replicate the logic if needed. Let's be consistent:
        total_size = end_frame - start_frame + 1
        all_data, adjusted_total_size = process_trajectory_slices_by_size(
            u, chain_output_dir, total_size, slice_size, chain_sele, start_frame
        )
        # Save the combined data
        save_data(all_data, chain_output_dir, trajectory_file, u, frames_used=adjusted_total_size)
    else:
        print("Error: You must specify either num_slices or slice_size.")
        raise RuntimeError("No slicing method specified.")

    # Collect CSV, RMSD, RMSF
    rmsx_csv = file_namer(chain_output_dir, trajectory_file, "csv", u=u, frames_used=adjusted_total_size)
    if verbose:
        print(f"RMSX CSV: {rmsx_csv}")
    rmsd_csv = calculate_rmsd(u, chain_output_dir, chain_sele=chain_sele, start_frame=start_frame,
                              end_frame=start_frame + adjusted_total_size - 1)
    if verbose and rmsd_csv:
        print(f"RMSD CSV: {rmsd_csv}")
    rmsf_csv = calculate_rmsf(u, chain_output_dir, chain_sele=chain_sele, start_frame=start_frame,
                              end_frame=start_frame + adjusted_total_size - 1)
    if verbose and rmsf_csv:
        print(f"RMSF CSV: {rmsf_csv}")

    # Update PDB B-factors with RMSX
    update_all_pdb_bfactors(rmsx_csv, silent=not verbose)

    # Conditionally make the plot now or skip for later
    if make_plot:
        if verbose:
            print("Generating plots...")
            print("This may take several minutes the first time it is run.")
        create_r_plot(
            rmsx_csv, rmsd_csv, rmsf_csv,
            rscript_executable, interpolate, triple, palette
        )
    else:
        if verbose:
            print("Skipping plot generation in run_rmsx() because make_plot=False.")

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
        sync_color_scale=False
):
    """
    Perform RMSX analysis for all chains in the topology file.

    Parameters
    ----------
    topology_file : str
        Path to the topology (e.g. .pdb, .gro) file.
    trajectory_file : str
        Path to the trajectory (e.g. .xtc, .dcd) file.
    output_dir : str, optional
        Path to the output directory where results will be written.
        If None, a default name based on the topology file is used.
    num_slices : int, optional
        Number of slices to split the trajectory into. If specified, slice_size is ignored.
    slice_size : int, optional
        Size (in frames) for each slice. If specified, num_slices is ignored.
    rscript_executable : str, optional
        Path to Rscript executable (default 'Rscript').
    verbose : bool, optional
        Enable detailed logging (default True).
    interpolate : bool, optional
        Enable interpolation for RMSX heatmaps (passed to geom_raster).
    triple : bool, optional
        If True, generate triple plots (RMSX, RMSD, RMSF) in addition to RMSX heatmaps.
    overwrite : bool, optional
        If True, overwrite the output directory if it exists. Otherwise prompt the user.
    palette : str, optional
        The color palette for the RMSX heatmap (e.g. 'viridis', 'plasma', etc.).
    start_frame : int, optional
        The first frame to include in the analysis. Default is 0.
    end_frame : int, optional
        The last frame to include (inclusive). Default None means all frames.
    sync_color_scale : bool, optional
        If True, defer plotting until all chains are analyzed, compute a global min/max
        of RMSX across all chains, and generate plots with that shared color scale.

    Returns
    -------
    combined_dir : str
        The final directory containing PDB files and/or plots. For multi-chain systems,
        this is typically <output_dir>/combined. For single-chain systems, it's the
        chain subfolder.

    Notes
    -----
    - If sync_color_scale=False, each chain is plotted and saved immediately based
      on its own min/max values.
    - If sync_color_scale=True, run_rmsx is called with make_plot=False for each chain,
      and then a global RMSX min/max is computed before plotting.
    - For single-chain cases, we point combined_dir to the chain subfolder so that
      subsequent tasks (e.g. run_flipbook) can find the correct PDB files.
    """
    import MDAnalysis as mda
    import numpy as np
    from pathlib import Path

    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)

    combined_output_dirs = []
    csv_paths = []

    for chain in chain_ids:
        print(f"\nAnalyzing Chain {chain}...")

        # If we want synced color scale, skip immediate plotting in run_rmsx
        chain_make_plot = not sync_color_scale

        run_rmsx(
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
            make_plot=chain_make_plot
        )

        chain_output_dir = os.path.join(output_dir, f"chain_{chain}_rmsx")
        combined_output_dirs.append(chain_output_dir)

        # Attempt to find the RMSX CSV in chain_output_dir
        possible_csv = list(Path(chain_output_dir).glob("rmsx_*.csv"))
        if possible_csv:
            csv_paths.append(str(possible_csv[0]))

    # Decide where we place the final PDB files (and final plots)
    if len(chain_ids) > 1:
        # Multi-chain: combine subfolders
        combined_dir = os.path.join(output_dir, "combined")
        print("\nCombining PDB files from all chains...")
        combine_pdb_files(combined_output_dirs, combined_dir)
        print("Combined RMSX analysis completed for all chains.")
    else:
        # Single chain: no need to combine multiple directories
        # Point combined_dir to the single chain subfolder
        single_chain_id = chain_ids[0]
        combined_dir = os.path.join(output_dir, f"chain_{single_chain_id}_rmsx")
        print(f"Single-chain analysis completed. Setting combined_dir to: {combined_dir}")

    # If sync_color_scale, compute global min/max and do the actual plotting now
    if sync_color_scale and csv_paths:
        print("\nComputing global RMSX min and max across all chains...")
        global_min, global_max = compute_global_rmsx_min_max(csv_paths)
        print(f"Global RMSX range = [{global_min:.3f}, {global_max:.3f}]")

        print("Generating plots for each chain with a fixed color scale...")
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
                max_val=global_max
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
        flipbook_min_bfactor=None,
        flipbook_max_bfactor=None,
        spacingFactor="1",
        start_frame=0,
        end_frame=None
):
    """
    Run RMSX analysis and generate a FlipBook visualization.
    """

    # The key change: set sync_color_scale=True
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
        sync_color_scale=True    # <-- Ensures consistent color scale across all chains
    )

    run_flipbook(
        directory=combined_dir,
        palette=palette,
        min_bfactor=flipbook_min_bfactor,
        max_bfactor=flipbook_max_bfactor,
        spacingFactor=spacingFactor
    )

    print("Full analysis including FlipBook visualization completed successfully.")

# The following sections are examples and tests. Uncomment and modify paths as needed to run.

# %%
# Example usage for single-chain RMSX analysis with uniform slice sizes
# traj_file_single = "/path/to/your/single_chain_trajectory.dcd"
# pdb_file_single = "/path/to/your/single_chain_protein.pdb"
# output_dir_single = "/path/to/your/output_directory_single_chain"

# run_rmsx_flipbook(
#     topology_file=pdb_file_single,
#     trajectory_file=traj_file_single,
#     output_dir=output_dir_single,
#     num_slices=16,             # Number of slices
#     slice_size=None,           # Not used since num_slices is specified
#     rscript_executable='Rscript',
#     verbose=True,
#     interpolate=False,
#     triple=True,
#     overwrite=True,
#     palette="magma",
#     spacingFactor="0.5",
#     start_frame=0,
#     end_frame=None  # Uses all frames up to the last, adjusted for uniform slices
# )

# %%
# Example usage for multi-chain RMSX analysis
# traj_file_multi = "/path/to/your/multi_chain_trajectory.xtc"
# pdb_file_multi = "/path/to/your/multi_chain_protein.pdb"
# output_dir_multi = "/path/to/your/output_directory_multi_chain"

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

# %%
# Additional comments and to-do items:
# should to test more md file types to make sure there aren't any issues with it telling the duration of the simulation
# could add hover over to show which time slice it came from, would need to add that info to flipbook somehow.
# allow you to say which frames you want from your simulation (allow you to cut it shorter)
# could allow it to be seen with RMSX per _ number of frames
# function that would auto find the files to display given the output_dir path - could be helpful for both the 3d plot but also for run_flipbook on existing dirs
# auto save and display image of
# log version of rmsx taking the log value of all the values in the rmsx df before ploting

# make flipbook have the option of running with a single chain... seems obvious in hindsight.
