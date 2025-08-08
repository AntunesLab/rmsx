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
import numpy as np  # Import numpy for unique operations
import pkg_resources # for managing the r scripts



#copied from moving_rmsx_utils_v10.py

def initialize_environment():
    """Print the Python executable path and current working directory."""
    print(sys.executable)
    print(os.getcwd())

def check_directory(output_dir):
    """Ask user if they want to overwrite an existing directory."""
    if os.path.exists(output_dir):
        response = input(f"The directory '{output_dir}' already exists. Do you want to overwrite it? (y/n): ")
        return response.strip().lower() == 'y'
    return True  # Proceed if the directory does not exist

def setup_directory(output_dir):
    """Set up or clear the output directory based on user preference."""
    if check_directory(output_dir):
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
        print(f"Process terminated by user. The directory '{output_dir}' will not be overwritten.")
        sys.exit(1)  # Exit the script if the user does not want to overwrite

def setup_universe(topology_file, trajectory_file):
    """Set up the MDAnalysis universe."""
    return mda.Universe(topology_file, trajectory_file)

def process_rmsx(psf_file, dcd_file, pdb_file, output_dir=None, slice_size=5, chain_sele=None):
    initialize_environment()
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        if chain_sele:
            output_dir = os.path.join(os.getcwd(), f"{base_name}_chain_{chain_sele}_rmsx")
        else:
            output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")
    # setup_directory(output_dir)
    u = setup_universe(psf_file, dcd_file)
    process_trajectory(u, pdb_file, dcd_file, output_dir, len(u.trajectory), slice_size, False, chain_sele)
    all_data = analyze_trajectory(output_dir, chain_sele)
    save_data(all_data, output_dir, dcd_file)

def process_trajectory(u, pdb_file, dcd_file, output_dir, total_size, slice_size, double_precision, chain_sele=None):
    """Process the trajectory by slicing it and saving slices, and write a PDB file for the first frame of each slice."""
    if total_size % slice_size != 0:
        print(f'LOSS OF FRAMES AT END! Total size {total_size} is not evenly divisible by slice size {slice_size}')

    u = mda.Universe(pdb_file, dcd_file)
    if u.trajectory.n_frames != total_size:
        raise ValueError(f'Unexpected number of frames: {u.trajectory.n_frames}, expected: {total_size}')
    print(f'The trajectory has {u.trajectory.n_frames} frames.')

    n_slices = total_size // slice_size
    for i in range(n_slices):
        start_frame = i * slice_size
        end_frame = (i + 1) * slice_size
        u.trajectory[start_frame]
        traj_path = os.path.join(output_dir, f'slice_{i + 1}.dcd')
        pdb_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')

        # Select protein atoms with chain selection
        selection_str = "protein"
        if chain_sele:
            selection_str += f" and segid {chain_sele}"
        protein = u.select_atoms(selection_str)

        # Write the first frame to PDB
        with mda.Writer(pdb_path, protein.n_atoms, multiframe=False) as pdb_writer:
            pdb_writer.write(protein)
        # Write the slice to DCD
        with mda.Writer(traj_path, protein.n_atoms) as dcd_writer:
            for ts in u.trajectory[start_frame:end_frame]:
                dcd_writer.write(protein)
        print(f'Slice {i + 1}: First frame saved to {pdb_path}, frames {start_frame} to {end_frame} saved to {traj_path}')

def analyze_trajectory(output_dir, chain_sele=None):
    """Analyze trajectory files by calculating RMSF and compiling results into a DataFrame."""
    all_data = pd.DataFrame()
    files_and_directories = os.listdir(output_dir)
    sorted_out = sorted(files_and_directories, key=lambda filename: [int(part) if part.isdigit() else part for part in re.split('(\d+)', filename)])

    for filename in sorted_out:
        if filename.endswith(".dcd"):
            pdb_filename = filename.replace('.dcd', '_first_frame.pdb')
            if pdb_filename in files_and_directories:
                dcd_filepath = os.path.join(output_dir, filename)
                pdb_filepath = os.path.join(output_dir, pdb_filename)
                print(f"Processing {dcd_filepath} with {pdb_filepath}")

                u = mda.Universe(pdb_filepath, dcd_filepath)
                selection_str = "protein and name CA"
                if chain_sele:
                    selection_str += f" and segid {chain_sele}"
                protein = u.select_atoms(selection_str)

                rmsf_analysis = RMSF(protein)
                rmsf_analysis.run()

                df = pd.DataFrame({filename: rmsf_analysis.rmsf}, index=[residue.resid for residue in protein.residues])
                all_data = pd.concat([all_data, df], axis=1) if not all_data.empty else df

    if not all_data.empty:
        u = mda.Universe(os.path.join(output_dir, sorted_out[0].replace('.dcd', '_first_frame.pdb')))
        selection_str = "protein and name CA"
        if chain_sele:
            selection_str += f" and segid {chain_sele}"
        protein = u.select_atoms(selection_str)
        all_data.insert(0, 'ChainID', [residue.atoms[0].segid for residue in protein.residues])
        all_data.insert(0, 'ResidueID', [residue.resid for residue in protein.residues])

    return all_data

def save_data(all_data, output_dir, dcd_file):
    """Save the compiled DataFrame to a CSV file."""
    output_filepath = file_namer(output_dir, dcd_file, "csv")
    # Save DataFrame to CSV
    all_data.to_csv(output_filepath, index=False)

def file_namer(output_dir, example_file, out_file_type="csv", prefix="rmsx"):
    simulation_length_fs = extract_simulation_length(example_file)
    sim_name = os.path.basename(example_file)
    sim_name = os.path.splitext(sim_name)[0]
    simulation_length_ns = simulation_length_fs / 1e6
    # Determine the number of decimal places based on the value
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

def extract_simulation_length(dcd_file):
    # Load the DCD file
    u = mda.Universe(dcd_file)
    timestep_fs = u.trajectory.dt * 1000.0  # Convert picoseconds to femtoseconds
    # Calculate total simulation length
    total_time_fs = timestep_fs * len(u.trajectory)
    return total_time_fs

def calculate_rmsd(pdb_file, dcd_file, output_dir, selection='protein and name CA', chain_sele=None):
    """
    Calculate the RMSD values for a given protein selection from a trajectory file against a reference structure.
    """
    u = mda.Universe(pdb_file, dcd_file)
    ref = mda.Universe(pdb_file)

    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)
    protein_ref = ref.select_atoms(selection)

    rmsd_analysis = rms.RMSD(protein, protein_ref)
    rmsd_analysis.run()

    columns = ['Frame', 'Time', 'RMSD']
    rmsd_df = pd.DataFrame(rmsd_analysis.rmsd, columns=columns)

    if output_dir:
        rmsd_output_filepath = os.path.join(output_dir, 'rmsd.csv')
        rmsd_df.to_csv(rmsd_output_filepath, index=False)
    return rmsd_output_filepath

def calculate_rmsf(pdb_file, dcd_file, output_dir=None, selection='protein and name CA', chain_sele=None):
    """
    Calculate the RMSF values for a given protein selection from a trajectory file.
    """
    u = mda.Universe(pdb_file, dcd_file)

    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)

    rmsf_analysis = RMSF(protein)
    rmsf_analysis.run()
    rmsf_values = rmsf_analysis.rmsf

    num_residues = len(protein.residues)
    print(f"Number of residues: {num_residues}")
    rmsf_list = rmsf_values.tolist()
    if len(rmsf_list) != num_residues:
        raise ValueError(f"Length of rmsf_list ({len(rmsf_list)}) does not match number of residues ({num_residues})")

    rmsf_whole_traj = pd.DataFrame({
        'ResidueID': [residue.resid for residue in protein.residues],
        'RMSF': rmsf_list
    })

    if output_dir:
        rmsf_output_filepath = os.path.join(output_dir, 'rmsf.csv')
        rmsf_whole_traj.to_csv(rmsf_output_filepath, index=False)
    return rmsf_output_filepath

def create_r_plot(rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable='Rscript', interpolate=False, triple=False):
    # r_script_path = os.path.expanduser(r_script_path)
    # Get the path to the R script within the package
    r_script_path = pkg_resources.resource_filename('rmsx', 'r_scripts/plot_rmsx.R')
    # Convert Python boolean to a string recognizable by R
    interpolate_str = 'TRUE' if interpolate is not None and interpolate else 'FALSE'
    triple_str = 'TRUE' if triple is not None and triple else 'FALSE'

    # Call the R script with subprocess
    result = subprocess.run([rscript_executable, r_script_path, rmsx_csv, rmsd_csv, rmsf_csv, interpolate_str, triple_str])

    if result.returncode == 0:
        print("R script executed successfully.")
    else:
        print("R script execution failed.")

    # Search for PNG files in the directory
    file_path = rmsx_csv
    directory = os.path.dirname(file_path)
    image_files = glob.glob(os.path.join(directory, '*.png'))

    # Display the first matched image
    if image_files:
        display(Image(filename=image_files[0]))
    else:
        print("No PNG files found in the specified directory.")

def run_rmsx(psf_file, dcd_file, pdb_file, output_dir=None, slice_size=5,
         rscript_executable='Rscript', verbose=True, interpolate=True, triple=False):
    # Initialize environment and set up output directory
    initialize_environment()
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")
        # Set up or clear the output directory based on user input
        # setup_directory(output_dir)

    # Load the PDB file to extract chain information
    u_pdb = mda.Universe(pdb_file)
    # Extract unique chain IDs (SEGIDs)
    chain_ids = np.unique(u_pdb.atoms.segids)
    # Get the number of residues in each chain
    chain_info = {}
    for chain in chain_ids:
        chain_atoms = u_pdb.select_atoms(f'segid {chain}')
        num_residues = len(chain_atoms.residues)
        chain_info[chain] = num_residues

    # Display available chains and their lengths
    print("Available chains and their lengths (in residues):")
    for chain, length in chain_info.items():
        print(f"Chain {chain}: {length} residues")

    # Prompt the user to select a chain
    chain_list = ", ".join([f"{chain} ({length} residues)" for chain, length in chain_info.items()])
    selected_chain = input(f"Please enter the chain ID you would like to analyze from the following options:\n{chain_list}\nChain ID: ").strip()

    # Check if the selected chain is valid
    if selected_chain not in chain_ids:
        print(f"Chain '{selected_chain}' is not available in the PDB file.")
        sys.exit(1)  # Exit the function or handle the error as desired

    # Proceed with analysis using the original PDB and DCD files, applying chain selection in analysis functions
    # Update output directory to include the chain ID
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_dir = os.path.join(output_dir, f"chain_{selected_chain}_rmsx")
    setup_directory(output_dir)

    # Run the analysis using the selected chain
    if verbose:
        print("Starting analysis...")
        process_rmsx(psf_file, dcd_file, pdb_file, output_dir, slice_size, chain_sele=selected_chain)
        rmsx_csv = file_namer(output_dir, dcd_file, "csv")
        print(f"RMSX CSV: {rmsx_csv}")
        rmsd_csv = calculate_rmsd(pdb_file, dcd_file, output_dir, chain_sele=selected_chain)
        print(f"RMSD CSV: {rmsd_csv}")
        rmsf_csv = calculate_rmsf(pdb_file, dcd_file, output_dir, chain_sele=selected_chain)
        print(f"RMSF CSV: {rmsf_csv}")
        print("Generating plots...")
        create_r_plot(rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple)
    else:
        with open(os.devnull, 'w') as f, redirect_stdout(f):
            process_rmsx(psf_file, dcd_file, pdb_file, output_dir, slice_size, chain_sele=selected_chain)
            rmsx_csv = file_namer(output_dir, dcd_file, "csv")
            rmsd_csv = calculate_rmsd(pdb_file, dcd_file, output_dir, chain_sele=selected_chain)
            rmsf_csv = calculate_rmsf(pdb_file, dcd_file, output_dir, chain_sele=selected_chain)
        with open(os.devnull, 'w') as f, redirect_stdout(f):
            create_r_plot(rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple)

