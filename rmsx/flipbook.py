#!/usr/bin/env python3

# prototype
# encapsulting flipbook to allow it to be imported into other programs more easily


"""
flipbook.py

A Python script to open all PDB files in a specified directory with ChimeraX,
sorted in natural numerical order based on the slice number in the filenames,
and apply additional ChimeraX commands such as dynamic coloring and tiling.

Usage:
    python3 flipbook.py /path/to/directory [--palette PALETTE] [--min_bfactor MIN] [--max_bfactor MAX]

Example:
    python3 flipbook.py /Users/finn/Desktop/RMSX_Demo_files_mac/gromacs_case_studies/new_protease_tfA/combined --palette viridis --min_bfactor 10 --max_bfactor 50

Integration with rmsx.py:
    You can import and use the run_flipbook function within rmsx.py to visualize results automatically.

    Example in rmsx.py:
        import flipbook

        # After analysis
        flipbook.run_flipbook(
            directory='/path/to/output_directory',
            palette='viridis',
            min_bfactor=10,
            max_bfactor=50
        )
"""

####


"""
flipbook.py

A Python script to open all PDB files in a specified directory with ChimeraX,
sorted in natural numerical order based on the slice number in the filenames,
and apply additional ChimeraX commands such as dynamic coloring and tiling.

"""

import os
import sys
import argparse
import subprocess
import re

# Available color palettes: "magma" (A), "inferno" (B), "plasma" (C), "viridis" (D), "cividis" (E), "rocket" (F), "mako" (G), "turbo" (H)

# Define color palettes with their respective hex codes
# The color hex codes were obtained from R's viridis_pal(option = ...)(12)
COLOR_PALETTES = {
    "magma": [  # From viridis_pal(option = "magma")(12)
        "#000004", "#120D32", "#331068", "#5A167E", "#7D2482",
        "#A3307E", "#C83E73", "#E95562", "#F97C5D", "#FEA873",
        "#FED395", "#FCFDBF"
    ],
    "inferno": [  # From viridis_pal(option = "inferno")(12)
        "#000004", "#140B35", "#3A0963", "#60136E", "#85216B",
        "#A92E5E", "#CB4149", "#E65D2F", "#F78311", "#FCAD12",
        "#F5DB4B", "#FCFFA4"
    ],
    "plasma": [  # From viridis_pal(option = "plasma")(12)
        "#0D0887", "#3E049C", "#6300A7", "#8707A6", "#A62098",
        "#C03A83", "#D5546E", "#E76F5A", "#F58C46", "#FDAD32",
        "#FCD225", "#F0F921"
    ],
    "viridis": [  # From viridis_pal(option = "viridis")(12)
        "#440154", "#482173", "#433E85", "#38598C", "#2D708E",
        "#25858E", "#1E9B8A", "#2BB07F", "#51C56A", "#85D54A",
        "#C2DF23", "#FDE725"
    ],
    "cividis": [  # From viridis_pal(option = "cividis")(12)
        "#00204D", "#00306F", "#2A406C", "#48526B", "#5E626E",
        "#727374", "#878479", "#9E9677", "#B6A971", "#D0BE67",
        "#EAD357", "#FFEA46"
    ],
    "rocket": [  # From viridis_pal(option = "rocket")(12)
        "#03051A", "#221331", "#451C47", "#6A1F56", "#921C5B",
        "#B91657", "#D92847", "#ED513E", "#F47C56", "#F6A47B",
        "#F7C9AA", "#FAEBDD"
    ],
    "mako": [  # From viridis_pal(option = "mako")(12)
        "#0B0405", "#231526", "#35264C", "#403A75", "#3D5296",
        "#366DA0", "#3487A6", "#35A1AB", "#43BBAD", "#6CD3AD",
        "#ADE3C0", "#DEF5E5"
    ],
    "turbo": [  # From viridis_pal(option = "turbo")(12)
        "#30123B", "#4454C4", "#4490FE", "#1FC8DE", "#29EFA2",
        "#7DFF56", "#C1F334", "#F1CA3A", "#FE922A", "#EA4F0D",
        "#BE2102", "#7A0403"
    ]
}


def create_color_mapping(palette_name, colors, min_bfactor, max_bfactor, num_models):
    """
    Creates a ChimeraX color byattribute command string based on the selected palette and B-factor range.

    Args:
        palette_name (str): The name of the color palette to use.
        colors (list): List of hex color codes.
        min_bfactor (float): The minimum B-factor value.
        max_bfactor (float): The maximum B-factor value.
        num_models (int): The number of models loaded.

    Returns:
        str: A formatted ChimeraX color byattribute command string.
    """
    num_colors = len(colors)
    if num_colors < 2:
        raise ValueError("At least two colors are required to create a gradient.")

    # Calculate the interval between color stops
    interval = (max_bfactor - min_bfactor) / (num_colors - 1)

    # Generate list of (bfactor_value, color) tuples
    color_stops = []
    for i, color in enumerate(colors):
        bfactor_value = round(min_bfactor + i * interval, 2)
        color_stops.append((bfactor_value, color))

    # Construct the palette mapping string with B-factor values first, then colors
    palette_mapping = ":".join([f"{bfactor},{color}" for bfactor, color in color_stops])

    # Construct the color byattribute command
    color_command = f"color byattribute a:bfactor #1-{num_models} target absc palette {palette_mapping}"

    return color_command


def extract_bfactor_range(pdb_file_paths):
    """
    Extracts the minimum and maximum B-factor values from the given PDB files.

    Args:
        pdb_file_paths (list): List of full file paths to PDB files.

    Returns:
        tuple: (min_bfactor, max_bfactor)
    """
    min_bfactor = float('inf')
    max_bfactor = float('-inf')

    for path in pdb_file_paths:
        try:
            with open(path, 'r') as pdb_file:
                for line in pdb_file:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            bfactor = float(line[60:66].strip())
                            if bfactor < min_bfactor:
                                min_bfactor = bfactor
                            if bfactor > max_bfactor:
                                max_bfactor = bfactor
                        except ValueError:
                            continue
        except Exception as e:
            print(f"Error reading file '{path}': {e}")
            continue

    # Handle cases where no B-factors were found
    if min_bfactor == float('inf') or max_bfactor == float('-inf'):
        min_bfactor, max_bfactor = 0.00, 1.00  # Default values

    return min_bfactor, max_bfactor


def find_pdb_files(directory, pattern=r'^slice_(\d+)_first_frame\.pdb$'):
    """
    Finds all PDB files in the specified directory matching the given regex pattern.

    Args:
        directory (str): The directory to search in.
        pattern (str): The regex pattern to match filenames.

    Returns:
        list: A list of matching filenames.
    """
    try:
        # List all entries in the directory
        entries = os.listdir(directory)
    except Exception as e:
        print(f"Error accessing directory '{directory}': {e}")
        sys.exit(1)

    # Compile the regex pattern
    regex = re.compile(pattern)

    # Filter files matching the regex
    pdb_files = [f for f in entries if os.path.isfile(os.path.join(directory, f)) and regex.match(f)]

    return pdb_files


def natural_sort_key(s):
    """
    Generates a key for natural sorting of strings containing numbers.

    Args:
        s (str): The string to generate a sort key for.

    Returns:
        list: A list of integers and lowercase strings for sorting.
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', s)]


def run_flipbook(directory, palette='viridis', min_bfactor=None, max_bfactor=None, spacingFactor=1,
                 extra_commands=None):
    """
    Executes the flipbook functionality to open PDB files in ChimeraX with specified settings.
    An optional extra_commands parameter allows you to append custom ChimeraX commands.

    Args:
        directory (str): Path to the directory containing PDB files.
        palette (str): Color palette to use for coloring the models.
        min_bfactor (float, optional): Minimum B-factor value. If None, auto-detect.
        max_bfactor (float, optional): Maximum B-factor value. If None, auto-detect.
        spacingFactor (int, optional): The spacing factor used in tiling.
        extra_commands (list or str, optional): Additional ChimeraX commands to run after the default commands.

    Returns:
        None
    """
    palette_name = palette
    provided_min_bfactor = min_bfactor
    provided_max_bfactor = max_bfactor

    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a valid directory.")
        sys.exit(1)

    pdb_files = find_pdb_files(directory)
    if not pdb_files:
        print(f"No files found matching pattern 'slice_<number>_first_frame.pdb' in directory '{directory}'.")
        sys.exit(1)

    pdb_files_sorted = sorted(pdb_files, key=natural_sort_key)
    pdb_file_paths = [os.path.join(directory, f) for f in pdb_files_sorted]
    num_models = len(pdb_file_paths)

    if provided_min_bfactor is None or provided_max_bfactor is None:
        min_bfactor, max_bfactor = extract_bfactor_range(pdb_file_paths)
        print(f"Detected B-factor range: {min_bfactor:.2f} - {max_bfactor:.2f}")
    else:
        min_bfactor = provided_min_bfactor
        max_bfactor = provided_max_bfactor
        print(f"Using provided B-factor range: {min_bfactor:.2f} - {max_bfactor:.2f}")

    colors = COLOR_PALETTES.get(palette_name)
    if not colors:
        print(f"Error: Palette '{palette_name}' is not defined.")
        sys.exit(1)

    try:
        color_command = create_color_mapping(palette_name, colors, min_bfactor, max_bfactor, num_models)
    except ValueError as ve:
        print(ve)
        sys.exit(1)

    columns = num_models
    axis_id = num_models + 1

    open_commands = " ; ".join([f"open '{path}'" for path in pdb_file_paths])

    default_commands = [
        "view",
        "define axis",
        f"view #{axis_id} zalign #{axis_id}",
        f"turn x 90 center #{axis_id}",
        "color byattribute bfactor",
        "worm bfactor",
        "lighting soft",
        "graphics silhouettes true",
        "set bgColor white",
        color_command,
        f"tile all columns {columns} spacingFactor {spacingFactor}",
        f"close #{axis_id}",
        f"save {directory}/rmsx_{palette}.png width 2000 height 1000 supersample 3 transparentBackground true"
    ]

    # Append custom extra commands, if provided.
    if extra_commands:
        if isinstance(extra_commands, str):
            extra_commands = [extra_commands]
        default_commands.extend(extra_commands)

    chimera_commands = f"{open_commands} ; " + " ; ".join(default_commands)
    cmd = ['chimerax', '--cmd', chimera_commands]

    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print("Error: 'chimerax' command not found. Please ensure ChimeraX is installed and added to your PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"ChimeraX exited with an error: {e}")
        sys.exit(e.returncode)
    except Exception as e:
        print(f"An unexpected error occurred while running ChimeraX: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Open PDB files in ChimeraX in numerical order and apply dynamic coloring and tiling.')
    parser.add_argument('directory', type=str, help='Path to the directory containing PDB files.')
    parser.add_argument('--palette', type=str, default='viridis', choices=COLOR_PALETTES.keys(),
                        help='Color palette to use for coloring the models.')
    parser.add_argument('--min_bfactor', type=float, default=None, help='Minimum B-factor value.')
    parser.add_argument('--max_bfactor', type=float, default=None, help='Maximum B-factor value.')
    # Allow extra commands to be specified on the command line; these will be appended.
    parser.add_argument('--extra-commands', type=str, nargs='+', default=[],
                        help='Extra ChimeraX commands to run after the default commands.')
    args = parser.parse_args()

    run_flipbook(
        directory=args.directory,
        palette=args.palette,
        min_bfactor=args.min_bfactor,
        max_bfactor=args.max_bfactor,
        extra_commands=args.extra_commands
    )


if __name__ == '__main__':
    main()
