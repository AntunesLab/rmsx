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
import pkg_resources  # for managing the R scripts
import plotly.graph_objects as go

def initialize_environment():
    """Print the Python executable path and current working directory."""
    print(sys.executable)
    print(os.getcwd())

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
        sys.exit(1)  # Exit the script if the user does not want to overwrite

def setup_universe(topology_file, trajectory_file):
    """Set up the MDAnalysis universe."""
    return mda.Universe(topology_file, trajectory_file)

def process_rmsx(topology_file, trajectory_file, output_dir=None, slice_size=5, chain_sele=None):
    initialize_environment()
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        if chain_sele:
            output_dir = os.path.join(
                os.getcwd(), f"{base_name}_chain_{chain_sele}_rmsx"
            )
        else:
            output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")
    u = setup_universe(topology_file, trajectory_file)
    process_trajectory(
        u,
        topology_file,
        trajectory_file,
        output_dir,
        len(u.trajectory),
        slice_size,
        False,
        chain_sele,
    )
    all_data = analyze_trajectory(output_dir, chain_sele)
    save_data(all_data, output_dir, trajectory_file)

def process_trajectory(
        u,
        topology_file,
        trajectory_file,
        output_dir,
        total_size,
        slice_size,
        double_precision,
        chain_sele=None,
):
    """Process the trajectory by slicing it and saving slices, and write a coordinate file for the first frame of each slice."""
    if total_size % slice_size != 0:
        print(
            f'LOSS OF FRAMES AT END! Total size {total_size} is not evenly divisible by slice size {slice_size}'
        )

    u = mda.Universe(topology_file, trajectory_file)
    if u.trajectory.n_frames != total_size:
        raise ValueError(
            f'Unexpected number of frames: {u.trajectory.n_frames}, expected: {total_size}'
        )
    print(f'The trajectory has {u.trajectory.n_frames} frames.')

    n_slices = total_size // slice_size
    for i in range(n_slices):
        start_frame = i * slice_size
        end_frame = (i + 1) * slice_size
        u.trajectory[start_frame]
        traj_path = os.path.join(output_dir, f'slice_{i + 1}.dcd')
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')

        # Select protein atoms with chain selection
        selection_str = "protein"
        if chain_sele:
            selection_str += f" and segid {chain_sele}"
        protein = u.select_atoms(selection_str)

        # Write the first frame to PDB
        with mda.Writer(coord_path, protein.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(protein)
        # Write the slice to DCD
        with mda.Writer(traj_path, protein.n_atoms) as dcd_writer:
            for ts in u.trajectory[start_frame:end_frame]:
                dcd_writer.write(protein)
        print(
            f'Slice {i + 1}: First frame saved to {coord_path}, frames {start_frame} to {end_frame} saved to {traj_path}'
        )

def analyze_trajectory(output_dir, chain_sele=None):
    """Analyze trajectory files by calculating RMSF and compiling results into a DataFrame."""
    all_data = pd.DataFrame()
    files_and_directories = os.listdir(output_dir)
    sorted_out = sorted(
        files_and_directories,
        key=lambda filename: [
            int(part) if part.isdigit() else part for part in re.split('(\d+)', filename)
        ],
    )

    for filename in sorted_out:
        if filename.endswith(".dcd"):
            coord_filename = filename.replace('.dcd', '_first_frame.pdb')
            if coord_filename in files_and_directories:
                dcd_filepath = os.path.join(output_dir, filename)
                coord_filepath = os.path.join(output_dir, coord_filename)
                print(f"Processing {dcd_filepath} with {coord_filepath}")

                u = mda.Universe(coord_filepath, dcd_filepath)
                selection_str = "protein and name CA"
                if chain_sele:
                    selection_str += f" and segid {chain_sele}"
                protein = u.select_atoms(selection_str)

                rmsf_analysis = RMSF(protein)
                rmsf_analysis.run()

                df = pd.DataFrame(
                    {filename: rmsf_analysis.rmsf},
                    index=[residue.resid for residue in protein.residues],
                )
                all_data = (
                    pd.concat([all_data, df], axis=1)
                    if not all_data.empty
                    else df
                )

    if not all_data.empty:
        u = mda.Universe(
            os.path.join(
                output_dir, sorted_out[0].replace('.dcd', '_first_frame.pdb')
            )
        )
        selection_str = "protein and name CA"
        if chain_sele:
            selection_str += f" and segid {chain_sele}"
        protein = u.select_atoms(selection_str)
        all_data.insert(
            0, 'ChainID', [residue.atoms[0].segid for residue in protein.residues]
        )
        all_data.insert(
            0, 'ResidueID', [residue.resid for residue in protein.residues]
        )

    return all_data

def save_data(all_data, output_dir, trajectory_file):
    """Save the compiled DataFrame to a CSV file."""
    output_filepath = file_namer(output_dir, trajectory_file, "csv")
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

def extract_simulation_length(trajectory_file):
    # Load the trajectory file
    u = mda.Universe(trajectory_file)
    timestep_fs = u.trajectory.dt * 1000.0  # Convert picoseconds to femtoseconds
    # Calculate total simulation length
    total_time_fs = timestep_fs * len(u.trajectory)
    return total_time_fs

def calculate_rmsd(
        topology_file,
        trajectory_file,
        output_dir,
        selection='protein and name CA',
        chain_sele=None,
):
    """
    Calculate the RMSD values for a given protein selection from a trajectory file against a reference structure.
    """
    u = mda.Universe(topology_file, trajectory_file)
    ref = mda.Universe(topology_file)

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

def calculate_rmsf(
        topology_file,
        trajectory_file,
        output_dir=None,
        selection='protein and name CA',
        chain_sele=None,
):
    """
    Calculate the RMSF values for a given protein selection from a trajectory file.
    """
    u = mda.Universe(topology_file, trajectory_file)

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
    return rmsf_output_filepath

def create_r_plot(
        rmsx_csv,
        rmsd_csv,
        rmsf_csv,
        rscript_executable='Rscript',
        interpolate=False,
        triple=False,
):
    # Get the path to the R script within the package
    r_script_path = pkg_resources.resource_filename(
        'rmsx', 'r_scripts/triple_plot_rmsx.R'
    )
    # Convert Python boolean to a string recognizable by R
    interpolate_str = 'TRUE' if interpolate else 'FALSE'
    triple_str = 'TRUE' if triple else 'FALSE'

    # Call the R script with subprocess
    result = subprocess.run(
        [
            rscript_executable,
            r_script_path,
            rmsx_csv,
            rmsd_csv,
            rmsf_csv,
            interpolate_str,
            triple_str,
        ]
    )

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

def update_pdb_bfactor(coord_file, rmsx_df):
    # Extract slice number from the file name
    coord_file_base_name = os.path.basename(coord_file)
    slice_number = int(
        coord_file_base_name.split('_')[1]
    )  # Assumes file name format is like 'slice_10_first_frame.pdb'
    # Create the column name for the current slice
    rmsx_column = f'slice_{slice_number}.dcd'
    # Open and read the coordinate file (PDB format)
    with open(coord_file, 'r') as pdb:
        pdb_lines = pdb.readlines()
    # Create a new list to store the modified PDB lines
    updated_pdb_lines = []
    for line in pdb_lines:
        # Check if the line is an ATOM or HETATM line
        if line.startswith(('ATOM', 'HETATM')):
            # Extract residue number (columns 23-26 in a PDB file)
            residue_number = int(line[22:26].strip())
            # Get the rmsx value from the DataFrame for this residue number using the slice column
            rmsx_value = rmsx_df.loc[
                rmsx_df['ResidueID'] == residue_number, rmsx_column
            ].values[0]
            # Update the B-factor (columns 61-66) with the rmsx value
            new_bfactor = f"{rmsx_value:6.2f}"  # Format as 6 characters, 2 decimal places
            updated_line = line[:60] + new_bfactor + line[66:]
            # Add the updated line to the list
            updated_pdb_lines.append(updated_line)
        else:
            # Non-ATOM/HETATM lines are copied as-is
            updated_pdb_lines.append(line)

    # Overwrite the original coordinate file with the updated lines
    with open(coord_file, 'w') as pdb:
        pdb.writelines(updated_pdb_lines)

    print(f"Original coordinate file {coord_file} has been updated with new B-factors.")

def load_coord_files(folder_path):
    # Function to extract numerical part of the file name
    def extract_number(filename):
        # Use regex to find numbers in the filename
        match = re.search(r'\d+', filename)
        return int(match.group()) if match else float('inf')  # Return infinity if no number is found
    # List all .pdb files in the folder
    coord_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    # Sort files based on the numerical value in the file name
    sorted_coord_files = sorted(coord_files, key=extract_number)
    # Return the full path of each file
    return [os.path.join(folder_path, f) for f in sorted_coord_files]

def update_all_pdb_bfactors(rmsx_csv):
    rmsx_df = pd.read_csv(rmsx_csv)
    coord_folder = os.path.dirname(rmsx_csv)
    coord_files = load_coord_files(coord_folder)
    for coord_file in coord_files:
        update_pdb_bfactor(coord_file, rmsx_df)

def plot_rmsx_surface(
        file_path: str,
        save_html: bool = False,
        html_filename: str = None,
        colorscale: str = "magma",
) -> None:
    """
    Generates an interactive 3D surface plot of RMSX data from a CSV file.

    Parameters:
    - file_path (str): Path to the RMSX data CSV file.
    - save_html (bool): If True, saves the plot as an HTML file. Default is False.
    - html_filename (str, optional):
        - If provided, specifies the filename for the saved HTML plot.
        - If None and `save_html` is True, the plot is saved as 'rmsx_surface_plot.html'
          in the same directory as the input CSV file.
    - colorscale (str): Plotly colorscale to use for the surface. Default is 'magma'.

    Returns:
    - None: Displays the interactive plot.
    """
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np
    import os

    # Resolve absolute path of the input file
    abs_file_path = os.path.abspath(file_path)
    input_dir = os.path.dirname(abs_file_path)

    # Load the data
    try:
        df = pd.read_csv(abs_file_path)
        df['ChainID'] = df['ChainID'].astype(str)  # Ensure ChainID is string
        print(f"Successfully loaded data from '{abs_file_path}'.")
    except FileNotFoundError:
        print(f"Error: The file '{abs_file_path}' was not found.")
        return
    except pd.errors.ParserError:
        print(
            f"Error: The file '{abs_file_path}' could not be parsed. Please check the file format."
        )
        return

    # Check if required columns exist
    required_columns = {'ResidueID', 'ChainID'}
    slice_columns = [col for col in df.columns if col.startswith('slice_')]
    if not required_columns.issubset(df.columns) or not slice_columns:
        print(
            "Error: The CSV file must contain 'ResidueID', 'ChainID', and 'slice_*.dcd' columns."
        )
        return

    # Get unique chains
    unique_chains = df['ChainID'].unique()
    if len(unique_chains) == 0:
        print("Error: No unique chains found in the data.")
        return

    # Prepare data for all chains
    chain_data = {}
    for chain in unique_chains:
        filtered = df[df['ChainID'] == chain]
        if filtered.empty:
            print(f"Warning: Chain '{chain}' has no data and will be skipped.")
            continue
        residue_ids = filtered['ResidueID'].values
        rmsx_matrix = filtered[slice_columns].values
        X, Y = np.meshgrid(np.arange(rmsx_matrix.shape[1]), residue_ids)
        chain_data[chain] = {'X': X, 'Y': Y, 'Z': rmsx_matrix}

    if not chain_data:
        print("Error: No valid chain data available for plotting.")
        return

    # Initialize the figure with the first chain
    initial_chain = unique_chains[0]
    if initial_chain not in chain_data:
        initial_chain = list(chain_data.keys())[0]

    # Convert initial_chain to string to avoid issues
    initial_chain_str = str(initial_chain)

    fig = go.Figure(
        data=[
            go.Surface(
                z=chain_data[initial_chain]['Z'],
                x=chain_data[initial_chain]['X'],
                y=chain_data[initial_chain]['Y'],
                colorscale=colorscale,  # Set colorscale to 'magma' by default
                name=f'Chain {initial_chain_str}',
                hovertemplate='Chain: %{name}<br>Slice: %{x}<br>Residue: %{y}<br>RMSX: %{z}<extra></extra>',
            )
        ]
    )

    # Define layout with dropdown menu for chain selection positioned at top-right
    fig.update_layout(
        title=f'3D RMSX Surface Plot for Chain {initial_chain_str}',
        scene=dict(
            xaxis_title='Slice',
            yaxis_title='Residue ID',
            zaxis_title='RMSX Value',
        ),
        autosize=True,
        width=1000,  # Increased width to provide more space for the dropdown
        height=800,
        margin=dict(
            l=65, r=250, b=65, t=90
        ),  # Increased right margin to accommodate dropdown
        updatemenus=[
            dict(
                buttons=[
                    dict(
                        args=[
                            {
                                'z': [chain_data[chain]['Z']],
                                'x': [chain_data[chain]['X']],
                                'y': [chain_data[chain]['Y']],
                                'title': f'3D RMSX Surface Plot for Chain {str(chain)}',
                                'name': f'Chain {str(chain)}',
                            }
                        ],
                        label=str(chain),  # Ensure label is string
                        method='update',
                    )
                    for chain in chain_data.keys()
                ],
                direction="down",
                showactive=True,
                x=1.15,  # Position further to the right
                xanchor="left",
                y=1,  # Align with the top of the plot
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",  # Semi-transparent background
                bordercolor="black",  # Border color
                borderwidth=1,  # Border width
            ),
        ],
    )

    # Display the plot
    fig.show()

    # Handle saving the plot as an HTML file
    if save_html:
        # Determine the HTML file path
        if html_filename is None:
            # Default filename in the same directory as the input file
            html_filename = "rmsx_surface_plot.html"

        # If html_filename is not an absolute path, place it in the input file's directory
        if not os.path.isabs(html_filename):
            abs_html_filename = os.path.join(input_dir, html_filename)
        else:
            abs_html_filename = html_filename

        try:
            fig.write_html(abs_html_filename)
            print(f"Plot successfully saved as '{abs_html_filename}'.")
        except Exception as e:
            print(f"Error: Failed to save HTML file. {e}")

def combine_pdb_files(chain_dirs, combined_dir):
    """
    Combines corresponding PDB files from multiple directories.

    Parameters:
        chain_dirs (list of str): A list of directories containing PDB files for different chains.
        combined_dir (str): The directory to store the combined PDB files.
    """
    # Create the combined output directory if it doesn't exist
    os.makedirs(combined_dir, exist_ok=True)

    # Assume all directories contain matching file names
    # Get a set of PDB filenames by checking the first directory
    pdb_files = [f for f in os.listdir(chain_dirs[0]) if f.endswith('.pdb')]

    for pdb_file in pdb_files:
        combined_content = []
        missing_files = False

        for chain_dir in chain_dirs:
            pdb_file_path = os.path.join(chain_dir, pdb_file)

            if os.path.exists(pdb_file_path):
                # Read the content of the current chain's PDB file
                with open(pdb_file_path, 'r') as file:
                    combined_content.extend(file.readlines())
            else:
                # If a file is missing in any chain directory, skip this file
                print(f"File {pdb_file} not found in {chain_dir}. Skipping.")
                missing_files = True
                break

        if not missing_files:
            # Write the combined content to the new file in the combined directory
            combined_pdb_file = os.path.join(combined_dir, pdb_file)
            with open(combined_pdb_file, 'w') as combined_file:
                combined_file.writelines(combined_content)

            print(f"Combined {pdb_file} from all chains into {combined_pdb_file}")

def find_and_combine_pdb_files(output_dir, combine_pdb_files_func):
    # Find directories ending with '_rmsx' and not starting with "combined"
    directories_with_dir = [
        os.path.join(output_dir, d)
        for d in os.listdir(output_dir)
        if d.endswith('_rmsx')
           and not d.startswith("combined")
           and os.path.isdir(os.path.join(output_dir, d))
    ]

    print("Directories to combine:", directories_with_dir)

    # Set up the combined directory path
    combined_dir = os.path.join(output_dir, "combined")

    # Call the combine_pdb_files function on the found directories
    combine_pdb_files_func(directories_with_dir, combined_dir)

def run_rmsx(
        topology_file,
        trajectory_file,
        output_dir=None,
        slice_size=5,
        rscript_executable='Rscript',
        verbose=True,
        interpolate=True,
        triple=False,
        chain_sele=None,
        overwrite=False,
):
    # Initialize environment and set up output directory
    initialize_environment()
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")

    # Load the topology file to extract chain information
    u_top = mda.Universe(topology_file)
    # Extract unique chain IDs (SEGIDs)
    chain_ids = np.unique(u_top.atoms.segids)
    # Get the number of residues in each chain
    chain_info = {}
    for chain in chain_ids:
        chain_atoms = u_top.select_atoms(f'segid {chain}')
        num_residues = len(chain_atoms.residues)
        chain_info[chain] = num_residues

    if chain_sele is None:
        # Display available chains and their lengths
        print("Available chains and their lengths (in residues):")
        for chain, length in chain_info.items():
            print(f"Chain {chain}: {length} residues")

        # Prompt the user to select a chain
        chain_list = ", ".join(
            [f"{chain} ({length} residues)" for chain, length in chain_info.items()]
        )
        selected_chain = input(
            f"Please enter the chain ID you would like to analyze from the following options:\n{chain_list}\nChain ID: "
        ).strip()

        # Check if the selected chain is valid
        if selected_chain not in chain_ids:
            print(f"Chain '{selected_chain}' is not available in the topology file.")
            sys.exit(1)  # Exit the function or handle the error as desired

        chain_sele = selected_chain  # Assign the selected chain
    else:
        # Check if the provided chain_sele is valid
        if chain_sele not in chain_ids:
            print(f"Chain '{chain_sele}' is not available in the topology file.")
            sys.exit(1)  # Exit the function or handle the error as desired

    # Update output directory to include the chain ID
    base_name = os.path.splitext(os.path.basename(topology_file))[0]
    output_dir = os.path.join(output_dir, f"chain_{chain_sele}_rmsx")
    setup_directory(output_dir, overwrite=overwrite)

    # Run the analysis using the selected chain
    if verbose:
        print("Starting analysis...")
        process_rmsx(topology_file, trajectory_file, output_dir, slice_size, chain_sele=chain_sele)
        rmsx_csv = file_namer(output_dir, trajectory_file, "csv")
        print(f"RMSX CSV: {rmsx_csv}")
        rmsd_csv = calculate_rmsd(
            topology_file, trajectory_file, output_dir, chain_sele=chain_sele
        )
        print(f"RMSD CSV: {rmsd_csv}")
        rmsf_csv = calculate_rmsf(
            topology_file, trajectory_file, output_dir, chain_sele=chain_sele
        )
        print(f"RMSF CSV: {rmsf_csv}")
        print("Generating plots...")
        print("This may take several minutes the first time it is run.")
        update_all_pdb_bfactors(rmsx_csv)
        create_r_plot(
            rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple
        )
    else:
        with open(os.devnull, 'w') as f, redirect_stdout(f):
            process_rmsx(
                topology_file, trajectory_file, output_dir, slice_size, chain_sele=chain_sele
            )
            rmsx_csv = file_namer(output_dir, trajectory_file, "csv")
            rmsd_csv = calculate_rmsd(
                topology_file, trajectory_file, output_dir, chain_sele=chain_sele
            )
            rmsf_csv = calculate_rmsf(
                topology_file, trajectory_file, output_dir, chain_sele=chain_sele
            )
            update_all_pdb_bfactors(rmsx_csv)
        with open(os.devnull, 'w') as f, redirect_stdout(f):
            create_r_plot(
                rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple
            )

def all_chain_rmsx(
        topology_file,
        trajectory_file,
        output_dir=None,
        slice_size=5,
        rscript_executable='Rscript',
        verbose=True,
        interpolate=True,
        triple=False,
        overwrite=False,
):
    u_top = mda.Universe(topology_file)
    # Extract unique chain IDs (SEGIDs)
    chain_ids = np.unique(u_top.atoms.segids)
    # Get the number of residues in each chain
    chain_info = {}
    for chain in chain_ids:
        chain_atoms = u_top.select_atoms(f'segid {chain}')
        num_residues = len(chain_atoms.residues)
        chain_info[chain] = num_residues

    for chain in chain_ids:
        run_rmsx(
            topology_file,
            trajectory_file,
            output_dir,
            slice_size,
            rscript_executable,
            verbose,
            interpolate,
            triple,
            chain_sele=chain,
            overwrite=overwrite,
        )

    print(f'Ran chains for {chain_ids}')

    find_and_combine_pdb_files(output_dir, combine_pdb_files)
    print("Combined all RMSX values across chains into one PDB per slice")
