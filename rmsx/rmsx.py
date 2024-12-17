
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
import pkg_resources
import plotly.graph_objects as go

from pathlib import Path
from .flipbook import run_flipbook  # Use a relative import


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
        raise RuntimeError("User chose not to overwrite the existing directory.")


def setup_universe(topology_file, trajectory_file):
    """Set up the MDAnalysis universe."""
    return mda.Universe(topology_file, trajectory_file)


def process_rmsx_by_slice_size(topology_file, trajectory_file, output_dir=None, slice_size=5, chain_sele=None):
    """
    Process RMSx by specifying the number of frames per slice.

    Parameters:
    - topology_file (str): Path to the topology file (e.g., PDB, PSF).
    - trajectory_file (str): Path to the trajectory file.
    - output_dir (str): Directory to save outputs. If None, defaults based on topology file name.
    - slice_size (int): Number of frames per slice.
    - chain_sele (str): Specific chain selection (e.g., 'A'). If None, selects all chains.
    """
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
    total_frames = len(u.trajectory)

    # Perform slicing and RMSF calculation
    all_data = process_trajectory_slices_by_size(u, output_dir, total_frames, slice_size, chain_sele)

    # Save the combined data
    save_data(all_data, output_dir, trajectory_file)


def process_trajectory_slices_by_size(u, output_dir, total_size, slice_size, chain_sele=None):
    """
    Slice the trajectory based on a fixed number of frames per slice.
    """
    if total_size % slice_size != 0:
        print(
            f'LOSS OF FRAMES AT END! Total size {total_size} is not evenly divisible by slice size {slice_size}'
        )

    n_slices = total_size // slice_size
    if n_slices == 0:
        raise ValueError("Slice size is larger than the total number of frames.")

    print(f'The trajectory has {u.trajectory.n_frames} frames.')
    print(f"Number of slices: {n_slices}")

    selection_str = "protein and name CA"
    if chain_sele:
        selection_str += f" and segid {chain_sele}"
    protein = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    for i in range(n_slices):
        start_frame = i * slice_size
        end_frame = (i + 1) * slice_size

        # Move to the start frame and write the first-frame PDB
        u.trajectory[start_frame]
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, protein.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(protein)
        print(f"First frame of slice {i + 1} written to {coord_path}")

        rmsf_calc = RMSF(protein)
        rmsf_calc.run(start=start_frame, stop=end_frame)

        df = pd.DataFrame(
            {f"slice_{i + 1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in protein.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        print(f'Slice {i + 1}: RMSF computed for frames {start_frame} to {end_frame - 1} (total {slice_size} frames)')

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(
            0, 'ChainID', [res.atoms[0].segid for res in protein.residues]
        )
        all_data.insert(
            0, 'ResidueID', [res.resid for res in protein.residues]
        )

    return all_data


def process_rmsx_by_num_slices(topology_file, trajectory_file, output_dir=None, num_slices=5, chain_sele=None):
    """
    Process RMSx by specifying the number of slices.
    """
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
    total_frames = len(u.trajectory)

    # Perform slicing and RMSF calculation
    all_data = process_trajectory_slices_by_num(u, output_dir, total_frames, num_slices, chain_sele)

    # Save the combined data
    save_data(all_data, output_dir, trajectory_file)


def process_trajectory_slices_by_num(u, output_dir, total_size, num_slices, chain_sele=None):
    """
    Slice the trajectory into a specific number of slices, distributing frames as evenly as possible.
    """
    if u.trajectory.n_frames != total_size:
        raise ValueError(
            f'Unexpected number of frames: {u.trajectory.n_frames}, expected: {total_size}'
        )
    print(f'The trajectory has {u.trajectory.n_frames} frames.')

    base_size = total_size // num_slices
    remainder = total_size % num_slices

    slice_sizes = [(base_size + 1) if i < remainder else base_size for i in range(num_slices)]

    print("Number of slices:", num_slices)
    for i, size in enumerate(slice_sizes, start=1):
        print(f"Slice {i}: {size} frames")

    selection_str = "protein and name CA"
    if chain_sele:
        selection_str += f" and segid {chain_sele}"
    protein = u.select_atoms(selection_str)

    all_data = pd.DataFrame()

    current_start = 0
    for i, size in enumerate(slice_sizes):
        start_frame = current_start
        end_frame = current_start + size  # end_frame is one past the last frame index for slicing
        current_start += size

        u.trajectory[start_frame]
        coord_path = os.path.join(output_dir, f'slice_{i + 1}_first_frame.pdb')
        with mda.Writer(coord_path, protein.n_atoms, multiframe=False) as coord_writer:
            coord_writer.write(protein)
        print(f"First frame of slice {i+1} written to {coord_path}")

        rmsf_calc = RMSF(protein)
        rmsf_calc.run(start=start_frame, stop=end_frame)

        df = pd.DataFrame(
            {f"slice_{i+1}.dcd": rmsf_calc.results.rmsf},
            index=[residue.resid for residue in protein.residues],
        )

        if all_data.empty:
            all_data = df
        else:
            all_data = pd.concat([all_data, df], axis=1)

        print(f'Slice {i + 1}: RMSF computed for frames {start_frame} to {end_frame - 1} (total {size} frames)')

    if not all_data.empty:
        # Add ResidueID and ChainID columns
        all_data.insert(
            0, 'ChainID', [res.atoms[0].segid for res in protein.residues]
        )
        all_data.insert(
            0, 'ResidueID', [res.resid for res in protein.residues]
        )

    return all_data


def analyze_trajectory(output_dir, chain_sele=None):
    # Not used in this approach
    return pd.DataFrame()


def save_data(all_data, output_dir, trajectory_file):
    output_filepath = file_namer(output_dir, trajectory_file, "csv")
    all_data.to_csv(output_filepath, index=False)


def file_namer(output_dir, example_file, out_file_type="csv", prefix="rmsx"):
    simulation_length_fs = extract_simulation_length(example_file)
    sim_name = os.path.basename(example_file)
    sim_name = os.path.splitext(sim_name)[0]
    simulation_length_ns = simulation_length_fs / 1e6
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
    u = mda.Universe(trajectory_file)
    timestep_fs = u.trajectory.dt * 1000.0  # ps to fs
    total_time_fs = timestep_fs * len(u.trajectory)
    return total_time_fs


def calculate_rmsd(
        topology_file,
        trajectory_file,
        output_dir,
        selection='protein and name CA',
        chain_sele=None,
):
    u = mda.Universe(topology_file, trajectory_file)
    ref = mda.Universe(topology_file)

    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)
    protein_ref = ref.select_atoms(selection)

    rmsd_analysis = rms.RMSD(protein, protein_ref)
    rmsd_analysis.run()

    columns = ['Frame', 'Time', 'RMSD']
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd, columns=columns)

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
    u = mda.Universe(topology_file, trajectory_file)

    if chain_sele:
        selection += f" and segid {chain_sele}"

    protein = u.select_atoms(selection)

    rmsf_analysis = RMSF(protein)
    rmsf_analysis.run()
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
    return rmsf_output_filepath


def create_r_plot(
        rmsx_csv,
        rmsd_csv,
        rmsf_csv,
        rscript_executable='Rscript',
        interpolate=False,
        triple=False,
        palette="plasma"
):
    interpolate_str = 'TRUE' if interpolate else 'FALSE'
    triple_str = 'TRUE' if triple else 'FALSE'

    try:
        try:
            current_dir = Path(__file__).parent.resolve()
        except NameError:
            current_dir = Path.cwd().resolve()

        r_script_path = current_dir / 'r_scripts' / 'plot_rmsx.R'

        if not r_script_path.is_file():
            print(f"Error: R script not found at {r_script_path}.")
            return

        print(f"Found R script at {r_script_path}.")

        result = subprocess.run(
            [
                rscript_executable,
                str(r_script_path),
                rmsx_csv,
                rmsd_csv,
                rmsf_csv,
                interpolate_str,
                triple_str,
                palette
            ],
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
            f"Error: Rscript executable not found: {rscript_executable}. Please ensure R is installed and 'Rscript' is in your PATH.")
        return
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return

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
    coord_file_base_name = os.path.basename(coord_file)
    slice_number = int(
        coord_file_base_name.split('_')[1]
    )
    rmsx_column = f'slice_{slice_number}.dcd'

    with open(coord_file, 'r') as pdb:
        pdb_lines = pdb.readlines()

    updated_pdb_lines = []
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            residue_number = int(line[22:26].strip())
            rmsx_value = rmsx_df.loc[
                rmsx_df['ResidueID'] == residue_number, rmsx_column
            ].values[0]
            new_bfactor = f"{rmsx_value:6.2f}"
            updated_line = line[:60] + new_bfactor + line[66:]
            updated_pdb_lines.append(updated_line)
        else:
            updated_pdb_lines.append(line)

    with open(coord_file, 'w') as pdb:
        pdb.writelines(updated_pdb_lines)
    if not silent:
        print(f"Original coordinate file {coord_file} has been updated with new B-factors.")


def load_coord_files(folder_path):
    def extract_number(filename):
        match = re.search(r'\d+', filename)
        return int(match.group()) if match else float('inf')

    coord_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    sorted_coord_files = sorted(coord_files, key=extract_number)
    return [os.path.join(folder_path, f) for f in sorted_coord_files]


def update_all_pdb_bfactors(rmsx_csv, silent):
    rmsx_df = pd.read_csv(rmsx_csv)
    coord_folder = os.path.dirname(rmsx_csv)
    coord_files = load_coord_files(coord_folder)
    for coord_file in coord_files:
        update_pdb_bfactor(coord_file, rmsx_df, silent)


def plot_rmsx_surface(
        file_path: str,
        save_html: bool = False,
        html_filename: str = None,
        colorscale: str = "magma",
) -> None:
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Plotly is not installed. Visualization features are disabled.")
        return

    import pandas as pd
    import numpy as np
    import os

    abs_file_path = os.path.abspath(file_path)
    input_dir = os.path.dirname(abs_file_path)

    try:
        df = pd.read_csv(abs_file_path)
        df['ChainID'] = df['ChainID'].astype(str)
        print(f"Successfully loaded data from '{abs_file_path}'.")
    except FileNotFoundError:
        print(f"Error: The file '{abs_file_path}' was not found.")
        return
    except pd.errors.ParserError:
        print(
            f"Error: The file '{abs_file_path}' could not be parsed. Please check the file format."
        )
        return

    required_columns = {'ResidueID', 'ChainID'}
    slice_columns = [col for col in df.columns if col.startswith('slice_')]
    if not required_columns.issubset(df.columns) or not slice_columns:
        print(
            "Error: The CSV file must contain 'ResidueID', 'ChainID', and 'slice_*.dcd' columns."
        )
        return

    unique_chains = df['ChainID'].unique()
    if len(unique_chains) == 0:
        print("Error: No unique chains found in the data.")
        return

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

    initial_chain = unique_chains[0]
    if initial_chain not in chain_data:
        initial_chain = list(chain_data.keys())[0]
    initial_chain_str = str(initial_chain)

    fig = go.Figure(
        data=[
            go.Surface(
                z=chain_data[initial_chain]['Z'],
                x=chain_data[initial_chain]['X'],
                y=chain_data[initial_chain]['Y'],
                colorscale=colorscale,
                name=f'Chain {initial_chain_str}',
                hovertemplate='Chain: %{name}<br>Slice: %{x}<br>Residue: %{y}<br>RMSX: %{z}<extra></extra>',
            )
        ]
    )

    fig.update_layout(
        title=f'3D RMSX Surface Plot for Chain {initial_chain_str}',
        scene=dict(
            xaxis_title='Slice',
            yaxis_title='Residue ID',
            zaxis_title='RMSX Value',
        ),
        autosize=True,
        width=1000,
        height=800,
        margin=dict(l=65, r=250, b=65, t=90),
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
                        label=str(chain),
                        method='update',
                    )
                    for chain in chain_data.keys()
                ],
                direction="down",
                showactive=True,
                x=1.15,
                xanchor="left",
                y=1,
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1,
            ),
        ],
    )

    fig.show()

    if save_html:
        if html_filename is None:
            html_filename = "rmsx_surface_plot.html"
        if not os.path.isabs(html_filename):
            abs_html_filename = os.path.join(input_dir, html_filename)
        else:
            abs_html_filename = html_filename
        try:
            fig.write_html(abs_html_filename)
            print(f"Plot successfully saved as '{abs_html_filename}'.")
        except Exception as e:
            print(f"Error: Failed to save HTML file. {e}")


def combine_pdb_files(chain_dirs, combined_dir, silent=False):
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
        palette="viridis"
):
    initialize_environment()
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(topology_file))[0]
        output_dir = os.path.join(os.getcwd(), f"{base_name}_rmsx")

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
    output_dir = os.path.join(output_dir, f"chain_{chain_sele}_rmsx")
    setup_directory(output_dir, overwrite=overwrite)

    if verbose:
        print("Starting analysis...")

        if num_slices is not None:
            # Use the new method: specify number of slices
            print(f"Using the new slicing method with num_slices={num_slices}")
            process_rmsx_by_num_slices(topology_file, trajectory_file, output_dir, num_slices, chain_sele=chain_sele)
        elif slice_size is not None:
            # Use the old method: specify slice size
            print(f"Using the old slicing method with slice_size={slice_size}")
            process_rmsx_by_slice_size(topology_file, trajectory_file, output_dir, slice_size, chain_sele=chain_sele)
        else:
            print("Error: You must specify either num_slices or slice_size.")
            raise RuntimeError("No slicing method specified.")

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

        update_all_pdb_bfactors(rmsx_csv, silent=False)

        print("Generating plots...")
        print("This may take several minutes the first time it is run.")
        create_r_plot(
            rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple, palette
        )
    else:
        with open(os.devnull, 'w') as f, redirect_stdout(f):
            if num_slices is not None:
                process_rmsx_by_num_slices(topology_file, trajectory_file, output_dir, num_slices,
                                           chain_sele=chain_sele)
            elif slice_size is not None:
                process_rmsx_by_slice_size(topology_file, trajectory_file, output_dir, slice_size,
                                           chain_sele=chain_sele)
            else:
                print("Error: You must specify either num_slices or slice_size.")
                raise RuntimeError("No slicing method specified.")

            rmsx_csv = file_namer(output_dir, trajectory_file, "csv")
            rmsd_csv = calculate_rmsd(
                topology_file, trajectory_file, output_dir, chain_sele=chain_sele
            )
            rmsf_csv = calculate_rmsf(
                topology_file, trajectory_file, output_dir, chain_sele=chain_sele
            )
            update_all_pdb_bfactors(rmsx_csv, silent=True)

        with open(os.devnull, 'w') as f, redirect_stdout(f):
            create_r_plot(
                rmsx_csv, rmsd_csv, rmsf_csv, rscript_executable, interpolate, triple, palette
            )


def all_chain_rmsx(
        topology_file,
        trajectory_file,
        output_dir=None,
        num_slices=None,  # Can specify either num_slices or slice_size
        slice_size=None,
        rscript_executable='Rscript',
        verbose=True,
        interpolate=True,
        triple=False,
        overwrite=False,
        palette='viridis'
):
    u_top = mda.Universe(topology_file)
    chain_ids = np.unique(u_top.atoms.segids)

    for chain in chain_ids:
        run_rmsx(
            topology_file,
            trajectory_file,
            output_dir,
            num_slices=num_slices,
            slice_size=slice_size,
            rscript_executable=rscript_executable,
            verbose=verbose,
            interpolate=interpolate,
            triple=triple,
            chain_sele=chain,
            overwrite=overwrite,
            palette=palette
        )

    print(f'Ran analysis for chains: {chain_ids}')

    if len(chain_ids) > 1:
        if verbose:
            find_and_combine_pdb_files(output_dir)
            print("Combined all RMSX values across chains into one PDB per slice")
        else:
            with open(os.devnull, 'w') as f, redirect_stdout(f):
                find_and_combine_pdb_files(output_dir)

        combined_dir = os.path.join(output_dir, "combined")
        return combined_dir
    else:
        directories_with_dir = [
            os.path.join(output_dir, d)
            for d in os.listdir(output_dir)
            if d.endswith('_rmsx')
               and not d.startswith("combined")
               and os.path.isdir(os.path.join(output_dir, d))
        ]
        print(f" Only one directory found, returning: {directories_with_dir[0]}")
        return directories_with_dir[0]


def run_rmsx_flipbook(
        topology_file,
        trajectory_file,
        output_dir=None,
        num_slices=9,  # Can specify either num_slices or slice_size, currently not implemented in the combined version?
        slice_size=None,
        rscript_executable='Rscript',
        verbose=True,
        interpolate=True,
        triple=False,
        overwrite=False,
        palette='viridis',
        flipbook_min_bfactor=None,
        flipbook_max_bfactor=None,
        spacingFactor="1"
):
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
        palette=palette
    )

    run_flipbook(
        directory=combined_dir,
        palette=palette,
        min_bfactor=flipbook_min_bfactor,
        max_bfactor=flipbook_max_bfactor,
        spacingFactor=spacingFactor
    )

    print("Full analysis including FlipBook visualization completed successfully.")

# should to test more md file types to make sure there aren't any issues with it telling the duration of the simulation
# could add hover over to show which time slice it came from
# allow you to say which frames you want from your simulation (allow you to cut it shorter)
# could allow it to be seen with RMSX per _ number of frames
# auto save picture in chimerax- last command could be visible so it could be copy pasted after it was correctly positioned.
# function that would auto find the files to display given the output_dir path - could be helpful for both the 3d plot but also for run_flipbook on existing dirs
# auto save and display image?
# log version of rmsx?
