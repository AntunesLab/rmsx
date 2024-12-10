import MDAnalysis as mda
import sys
import os

def detect_file_format(file_path):
    """
    Detects the file format based on the file extension.
    Returns the file format string.
    """
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()
    # Mapping of extensions to MDAnalysis format names
    format_map = {
        '.dcd': 'DCD',
        '.xtc': 'XTC',
        '.trr': 'TRR',
        '.nc': 'NCDF',
        '.ncdf': 'NCDF',
        '.pdb': 'PDB',
        '.xyz': 'XYZ',
        '.lammpstrj': 'LAMMPSDUMP',
        '.gsd': 'GSD',
        '.trj': 'TRJ',
        '.mdcrd': 'CRD',
        '.rst': 'RESTART',
        '.gro': 'GRO',
        '.h5md': 'H5MD'
    }
    return format_map.get(ext, None)

def guess_simulation_engine(file_format):
    """
    Guess the simulation engine based on the file format.
    This is a heuristic and may not be 100% accurate.
    """
    # Common associations:
    # DCD: often from CHARMM/NAMD (but can also come from other sources)
    # XTC/TRR/GRO: GROMACS
    # LAMMPSDUMP: LAMMPS
    # GSD: HOOMD-blue
    # NCDF: AMBER
    engine_map = {
        'DCD': 'CHARMM/NAMD',
        'XTC': 'GROMACS',
        'TRR': 'GROMACS',
        'GRO': 'GROMACS',
        'NCDF': 'AMBER',
        'LAMMPSDUMP': 'LAMMPS',
        'GSD': 'HOOMD-blue'
    }
    return engine_map.get(file_format, 'Unknown engine')

def get_time_unit(file_format):
    """
    Returns the time unit for the given file format.
    If the time unit is unknown or variable, returns None.
    """
    # Mapping of file formats to time units and conversion factors to picoseconds
    # If conversion factor is None, time is assumed to be in picoseconds
    time_unit_map = {
        'DCD': ('AKMA', 0.04888821),         # CHARMM/NAMD DCD files use AKMA units
        'LAMMPSDUMP': (None, None),          # LAMMPS units are user-defined
        'GSD': ('reduced', None),            # HOOMD-blue uses reduced units
        'XTC': ('picoseconds', None),        # GROMACS typically uses picoseconds
        'TRR': ('picoseconds', None),        # GROMACS typically uses picoseconds
        'NCDF': ('picoseconds', None),       # AMBER NetCDF files usually use picoseconds
        'H5MD': (None, None),                # H5MD units can vary
    }
    return time_unit_map.get(file_format, (None, None))

def convert_time_to_ps(time, unit):
    """
    Converts time to picoseconds based on the given unit.
    Returns time in picoseconds.
    """
    # Define conversion factors to picoseconds
    conversion_factors = {
        'AKMA': 0.04888821,
        'femtoseconds': 0.001,
        'picoseconds': 1.0,
        'nanoseconds': 1000.0,
        'reduced': None,  # User needs to provide conversion factor
    }
    factor = conversion_factors.get(unit)
    if factor is None:
        return None  # Unknown unit, cannot convert
    return time * factor

def calculate_simulation_details(u, file_format):
    """
    Calculates simulation details:
    - Simulation engine (based on file_format)
    - Timestep per frame (in ps)
    - Total simulation duration (in ns)
    """
    sim_engine = guess_simulation_engine(file_format)
    time_unit, conversion_factor = get_time_unit(file_format)

    # Prompt user if time_unit is unknown
    if time_unit is None:
        time_unit = input(f"Enter the time unit for the {file_format} trajectory (e.g., picoseconds, femtoseconds): ").strip()
        if time_unit not in ['femtoseconds', 'picoseconds', 'nanoseconds']:
            print("Unsupported or unknown time unit. Exiting.")
            sys.exit(1)

    times = [ts.time for ts in u.trajectory]

    # Check if time info is available
    if all(t == 0 for t in times):
        # Time not in trajectory, try dt or ask user
        timestep = u.trajectory.dt
        if timestep is None or timestep == 0:
            try:
                timestep = float(input(f"Enter the timestep between frames (in {time_unit}): "))
            except ValueError:
                print("Invalid input. Exiting.")
                sys.exit(1)
    else:
        # Calculate timestep from times array
        if len(times) >= 2:
            timestep = times[1] - times[0]
        else:
            # Only one frame? Unclear timestep, ask user
            try:
                timestep = float(input(f"Only one frame detected. Enter the timestep (in {time_unit}): "))
            except ValueError:
                print("Invalid input. Exiting.")
                sys.exit(1)

    timestep_ps = convert_time_to_ps(timestep, time_unit)
    if timestep_ps is None:
        # Ask user for conversion factor
        try:
            user_factor = float(input(f"Enter the conversion factor from {time_unit} to picoseconds: "))
            timestep_ps = timestep * user_factor
        except ValueError:
            print("Invalid conversion factor. Exiting.")
            sys.exit(1)

    num_frames = len(u.trajectory)
    total_duration_ps = timestep_ps * (num_frames - 1)
    total_duration_ns = total_duration_ps / 1000.0

    return sim_engine, timestep_ps, total_duration_ns, num_frames

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Extract simulation details from an MD trajectory.")
    parser.add_argument('-t', '--topology', help='Topology file (e.g., PSF, PDB)', required=False)
    parser.add_argument('-traj', '--trajectory', help='Trajectory file', required=True)
    parser.add_argument('-manual', '--manual_duration', help='Manually specify the simulation duration in nanoseconds', type=float)
    args = parser.parse_args()

    topology_file = args.topology
    trajectory_file = args.trajectory
    manual_duration = args.manual_duration

    # Detect file format
    file_format = detect_file_format(trajectory_file)
    if file_format is None:
        print("Unsupported or unknown trajectory file format.")
        sys.exit(1)
    else:
        print(f"Detected trajectory format: {file_format}")

    # Load the Universe
    try:
        if topology_file:
            u = mda.Universe(topology_file, trajectory_file, format=file_format)
        else:
            u = mda.Universe(trajectory_file, format=file_format)
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        sys.exit(1)

    if manual_duration is not None:
        total_duration_ns = manual_duration
        print(f"Using manually specified simulation duration: {total_duration_ns} ns")
        # Cannot infer timestep or engine easily if we skip calculation
        # If needed, we could re-run calculation anyway or just print partial info
        sim_engine = guess_simulation_engine(file_format)
        print(f"Possible Simulation Engine: {sim_engine}")
    else:
        sim_engine, timestep_ps, total_duration_ns, num_frames = calculate_simulation_details(u, file_format)
        print(f"Simulation Engine Guess: {sim_engine}")
        print(f"Timestep per frame: {timestep_ps:.4f} ps")
        print(f"Total number of frames: {num_frames}")
        print(f"Total simulation duration: {total_duration_ns:.3f} ns")
        print(f"Note: One frame corresponds to {timestep_ps:.4f} ps of simulation time.")

# if __name__ == "__main__":
#     main()
