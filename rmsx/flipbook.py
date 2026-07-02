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
import platform
import shutil
import csv
from pathlib import Path
from collections import defaultdict

try:
    import pty  # type: ignore[attr-defined]
except Exception:  # pragma: no cover (not available on Windows)
    pty = None

import threading

from rmsx.vmd_scripts.vmd_finder import find_vmd_executable

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
        "#0B0405", "#231526", "#35264C", "#403A75", "#3D526D",
        "#366DA0", "#3487A6", "#35A1AB", "#43BBAD", "#6CD3AD",
        "#ADE3C0", "#DEF5E5"
    ],
    "turbo": [  # From viridis_pal(option = "turbo")(12)
        "#30123B", "#4454C4", "#4490FE", "#1FC8DE", "#29EFA2",
        "#7DFF56", "#C1F334", "#F1CA3A", "#FE922A", "#EA4F0D",
        "#BE2102", "#7A0403"
    ]
}

MASK_METADATA_FILENAME = "masked_residues.csv"
MASK_TRANSPARENCY_PERCENT = 70
VMD_MASK_OPACITY = 0.30

# --------------------------- ChimeraX discovery ---------------------------

def _which_chimerax(explicit=None):
    """
    Locate a ChimeraX executable across platforms.
    Resolution order:
      1) explicit path argument
      2) CHIMERAX_EXECUTABLE env var
      3) PATH: 'ChimeraX'/'chimerax'/'ChimeraX.exe'
      4) OS-specific well-known locations
         - macOS: pick newest ChimeraX*.app or UCSF ChimeraX*.app in /Applications or ~/Applications
         - Linux: /usr/local/bin/chimerax, /usr/bin/chimerax, /opt/ChimeraX/bin/chimerax
         - Windows: %ProgramFiles%/ChimeraX*/bin/ChimeraX.exe (incl. versioned folders)
    Prints the resolved path when found.
    """
    def add(label, p, sources, candidates):
        sources.append((label, p))
        candidates.append(p)

    sources, candidates = [], []

    # 1) explicit
    if explicit:
        add("arg --chimerax", Path(explicit), sources, candidates)

    # 2) env
    env = os.environ.get("CHIMERAX_EXECUTABLE")
    if env:
        add("env CHIMERAX_EXECUTABLE", Path(env), sources, candidates)

    # 3) PATH
    for name in ("ChimeraX", "chimerax", "ChimeraX.exe"):
        hit = shutil.which(name)
        if hit:
            add(f"PATH:{name}", Path(hit), sources, candidates)

    # 4) OS-specific
    system = platform.system()

    if system == "Darwin":
        import re as _re

        def find_mac_apps(base: Path):
            out = []
            if base.exists():
                for pat in ("ChimeraX*.app", "UCSF ChimeraX*.app"):
                    for app in base.glob(pat):
                        m = _re.search(r"(\d+(?:\.\d+){0,3})", app.name)
                        vtuple = tuple(int(x) for x in m.group(1).split(".")) if m else (0,)
                        is_daily = "daily" in app.name.lower()
                        exe = app / "Contents" / "MacOS" / "ChimeraX"
                        out.append((vtuple, is_daily, app.name, exe))
            return out

        versions = []
        versions += find_mac_apps(Path("/Applications"))
        versions += find_mac_apps(Path.home() / "Applications")

        if versions:
            versions.sort(key=lambda t: (t[0], not t[1]), reverse=True)
            print("[flipbook] macOS bundle candidates (sorted):")
            for v, is_daily, name, exe in versions:
                vv = ".".join(map(str, v)) if v != (0,) else "0"
                print(f"   - {name} -> {exe}  (ver={vv}, daily={is_daily})")
            _, _, appname, execp = versions[0]
            add(f"macOS app bundle ({appname})", execp, sources, candidates)
        else:
            print("[flipbook] macOS: no ChimeraX app bundles found in /Applications or ~/Applications")

    elif system == "Linux":
        for guess in ("/usr/local/bin/chimerax", "/usr/bin/chimerax", "/opt/ChimeraX/bin/chimerax"):
            add(f"linux default:{guess}", Path(guess), sources, candidates)

    elif system == "Windows":
        pf  = os.environ.get("ProgramFiles", r"C:\Program Files")
        pfx = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
        guesses = [
            Path(pf)  / "ChimeraX" / "bin" / "ChimeraX.exe",
            Path(pfx) / "ChimeraX" / "bin" / "ChimeraX.exe",
        ]
        for base in (Path(pf), Path(pfx)):
            for d in base.glob("ChimeraX*"):
                guesses.append(d / "bin" / "ChimeraX.exe")
        for g in guesses:
            add(f"windows default:{g}", g, sources, candidates)

    print("[flipbook] Checking ChimeraX candidates (in order):")
    for label, path in sources:
        print(f"   - {label}: {path}")
        if path.exists() and os.access(path, os.X_OK):
            print(f"[flipbook] Using ChimeraX installation: {path}")
            return str(path)

    raise FileNotFoundError(
        "Could not find ChimeraX. "
        "Set CHIMERAX_EXECUTABLE, pass an explicit path, or add it to PATH."
    )

# --------------------------- coloring helpers -----------------------------

def create_color_mapping(palette_name, colors, min_bfactor, max_bfactor, num_models):
    """
    Creates a ChimeraX color byattribute command string based on the selected palette and B-factor range.
    """
    num_colors = len(colors)
    if num_colors < 2:
        raise ValueError("At least two colors are required to create a gradient.")

    interval = (max_bfactor - min_bfactor) / (num_colors - 1)

    color_stops = []
    for i, color in enumerate(colors):
        bfactor_value = round(min_bfactor + i * interval, 2)
        color_stops.append((bfactor_value, color))

    palette_mapping = ":".join([f"{bfactor},{color}" for bfactor, color in color_stops])

    return f"color byattribute a:bfactor #1-{num_models} target absc palette {palette_mapping}"


def extract_bfactor_range(pdb_file_paths):
    """
    Extracts the minimum and maximum B-factor values from the given PDB files.
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

    if min_bfactor == float('inf') or max_bfactor == float('-inf'):
        min_bfactor, max_bfactor = 0.00, 1.00

    return min_bfactor, max_bfactor


def find_pdb_files(directory, pattern=r'^slice_(\d+)_first_frame\.pdb$'):
    """
    Finds all PDB files in the specified directory matching the given regex pattern.
    """
    try:
        entries = os.listdir(directory)
    except Exception as e:
        raise OSError(f"Error accessing directory '{directory}': {e}") from e

    regex = re.compile(pattern)
    return [f for f in entries if os.path.isfile(os.path.join(directory, f)) and regex.match(f)]


def natural_sort_key(s):
    """
    Natural sorting for slice numbers.
    """
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def _parse_mask_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _load_masked_residue_map(directory):
    mask_path = Path(directory) / MASK_METADATA_FILENAME
    if not mask_path.is_file():
        return {}

    masked_map = defaultdict(set)
    with mask_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if not _parse_mask_bool(row.get("Masked", "")):
                continue
            chain_id = str(row.get("ChainID", "")).strip()
            residue_text = str(row.get("ResidueID", "")).strip()
            if not residue_text:
                continue
            try:
                residue_id = int(float(residue_text))
            except ValueError:
                continue
            masked_map[chain_id].add(residue_id)
    return {chain_id: sorted(residues) for chain_id, residues in masked_map.items()}


def _build_pdb_chain_lookup(pdb_file_path):
    chain_lookup = {}
    with open(pdb_file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain_id = line[21].strip()
            segid = line[72:76].strip()
            if segid and segid not in chain_lookup:
                chain_lookup[segid] = chain_id or segid
            if chain_id and chain_id not in chain_lookup:
                chain_lookup.setdefault(chain_id, chain_id)
    return chain_lookup


def _group_residue_ranges(residue_ids):
    if not residue_ids:
        return []
    sorted_ids = sorted(set(int(residue_id) for residue_id in residue_ids))
    ranges = []
    start_id = sorted_ids[0]
    end_id = sorted_ids[0]
    for residue_id in sorted_ids[1:]:
        if residue_id == end_id + 1:
            end_id = residue_id
            continue
        ranges.append((start_id, end_id))
        start_id = residue_id
        end_id = residue_id
    ranges.append((start_id, end_id))
    return ranges


def _format_residue_ranges(residue_ids):
    range_strings = []
    for start_id, end_id in _group_residue_ranges(residue_ids):
        if start_id == end_id:
            range_strings.append(str(start_id))
        else:
            range_strings.append(f"{start_id}-{end_id}")
    return ",".join(range_strings)


def _normalize_chimerax_chain_id(chain_id):
    stripped = str(chain_id).strip()
    if not stripped:
        return ""
    if len(stripped) == 1 and stripped.isalpha():
        return stripped
    if len(stripped) == 1 and stripped.isdigit():
        return "X"
    return "X"


def _build_chimerax_mask_commands(directory, pdb_file_paths):
    masked_residues = _load_masked_residue_map(directory)
    if not masked_residues:
        return []

    chain_lookup = _build_pdb_chain_lookup(pdb_file_paths[0])
    transparency_commands = []
    for original_chain_id, residue_ids in sorted(masked_residues.items()):
        if not residue_ids:
            continue
        pdb_chain_id = chain_lookup.get(original_chain_id, original_chain_id)
        chimerax_chain_id = _normalize_chimerax_chain_id(pdb_chain_id)
        residue_spec = _format_residue_ranges(residue_ids)
        if not residue_spec:
            continue
        if chimerax_chain_id:
            target_spec = f"/{chimerax_chain_id}:{residue_spec}"
        else:
            target_spec = f":{residue_spec}"
        transparency_commands.append(
            f"transparency {target_spec} {MASK_TRANSPARENCY_PERCENT} target r"
        )
    return transparency_commands


# -------------------------------------------------------------------------
# ------------------------- Flipbook MAIN FUNCTION -------------------------
# -------------------------------------------------------------------------

def run_flipbook(directory, palette='viridis', min_bfactor=None, max_bfactor=None,
                 spacingFactor=None, extra_commands=None, viewer='chimerax',
                 molstar_output=None, molstar_manifest_output=None,
                 molstar_asset_mode='cdn', molstar_height=720,
                 molstar_camera_mode='orthographic'):
    """
    Executes the flipbook functionality to open PDB files in the chosen viewer.
    """
    palette_name = palette
    provided_min_bfactor = min_bfactor
    provided_max_bfactor = max_bfactor

    if not os.path.isdir(directory):
        raise NotADirectoryError(f"'{directory}' is not a valid directory.")

    pdb_files = find_pdb_files(directory)
    if not pdb_files:
        raise FileNotFoundError(
            f"No files found matching pattern 'slice_<number>_first_frame.pdb' in directory '{directory}'."
        )

    pdb_files_sorted = sorted(pdb_files, key=natural_sort_key)
    pdb_file_paths = [os.path.join(directory, f) for f in pdb_files_sorted]
    num_models = len(pdb_file_paths)
    mask_commands = _build_chimerax_mask_commands(directory, pdb_file_paths)
    mask_metadata_path = Path(directory) / MASK_METADATA_FILENAME
    has_mask_metadata = mask_metadata_path.is_file() and bool(_load_masked_residue_map(directory))

    if provided_min_bfactor is None or provided_max_bfactor is None:
        min_bfactor, max_bfactor = extract_bfactor_range(pdb_file_paths)
        print(f"Detected B-factor range: {min_bfactor:.2f} - {max_bfactor:.2f}")
    else:
        min_bfactor = provided_min_bfactor
        max_bfactor = provided_max_bfactor
        print(f"Using provided B-factor range: {min_bfactor:.2f} - {max_bfactor:.2f}")

    colors = COLOR_PALETTES.get(palette_name)
    if not colors:
        raise ValueError(f"Palette '{palette_name}' is not defined.")

    try:
        color_command = create_color_mapping(palette_name, colors, min_bfactor, max_bfactor, num_models)
    except ValueError as ve:
        raise ValueError(str(ve)) from ve

    columns = num_models
    axis_id = num_models + 1

    open_commands = " ; ".join([f"open '{path}'" for path in pdb_file_paths])

    chimerax_spacing_factor = 1 if spacingFactor is None else spacingFactor

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
        f"tile all columns {columns} spacingFactor {chimerax_spacing_factor}",
        f"close #{axis_id}",
        f"save {directory}/rmsx_{palette}.png width 2000 height 1000 supersample 3 transparentBackground true"
    ]

    if mask_commands:
        insertion_index = default_commands.index("lighting soft")
        default_commands[insertion_index:insertion_index] = mask_commands

    if extra_commands:
        if isinstance(extra_commands, str):
            extra_commands = [extra_commands]
        default_commands.extend(extra_commands)

    #
    # =====================================================================
    # --------------------------- VMD BACKEND -----------------------------
    # =====================================================================
    #
    viewer_name = viewer.lower()

    if viewer_name == "molstar":
        from .molstar_viewer import write_molstar_flipbook

        result = write_molstar_flipbook(
            directory=directory,
            palette=palette_name,
            min_bfactor=min_bfactor,
            max_bfactor=max_bfactor,
            spacing_factor=spacingFactor,
            camera_mode=molstar_camera_mode,
            output_html=molstar_output,
            output_manifest=molstar_manifest_output,
            asset_mode=molstar_asset_mode,
            iframe_height=molstar_height,
            title="RMSX Molstar Flipbook",
            mask_filename=MASK_METADATA_FILENAME,
            color_palettes=COLOR_PALETTES,
        )
        print(f"[flipbook] Molstar HTML written to: {result.html_path}")
        print(f"[flipbook] Molstar manifest written to: {result.manifest_path}")
        return result

    if viewer_name == "vmd":
        print("[flipbook] Using VMD backend (PTY Mode)...")

        # 1. Define the PTY reader thread function
        #    This runs in the background to read VMD's output
        def pty_reader(fd):
            """Reads output from the PTY's master file descriptor."""
            try:
                while True:
                    data_b = os.read(fd, 1024)
                    if not data_b:
                        break  # PTY closed
                    # Optional: uncomment below to see VMD's live output
                    # print(f"[VMD PTY]: {data_b.decode('utf-8', errors='ignore').rstrip()}")
            except Exception:
                pass # This will error when the process closes, which is normal
            finally:
                os.close(fd)

        vmd_loader = os.path.join(
            os.path.dirname(__file__),
            "vmd_scripts",
            "wait_to_load.tcl"
        )
        # "grid_color_scale_centered_xaxis_hotkeys.tcl"
        if not os.path.exists(vmd_loader):
            raise FileNotFoundError(f"VMD loader script not found at: {vmd_loader}")

        try:
            vmd_exec = find_vmd_executable()

            # Build the simple, original command
            cmd = [
                vmd_exec,
                "-dispdev", "win",
                "-e", vmd_loader,
                "-args",
                *pdb_file_paths,
                palette_name]

            print(f"[flipbook] Launching VMD in PTY:\n  {' '.join(cmd)}")

            vmd_cwd = os.path.dirname(vmd_loader)
            vmd_main = os.path.join(vmd_cwd, "grid_color_scale_centered_xaxis_hotkeys.tcl")
            vmd_env = os.environ.copy()
            vmd_env["RMSX_VMD_MAIN"] = vmd_main
            vmd_env.setdefault("RMSX_VMD_THICKFIELD", "user")
            if has_mask_metadata:
                vmd_env["RMSX_VMD_MASK_FILE"] = str(mask_metadata_path)
                vmd_env["RMSX_VMD_MASK_OPACITY"] = str(VMD_MASK_OPACITY)
            if pty is None:
                # Windows (and some restricted environments) do not support PTYs.
                # Fall back to a detached process without attaching to a pseudo-terminal.
                subprocess.Popen(cmd, start_new_session=True, cwd=vmd_cwd, env=vmd_env)
            else:
                # 2. Create the pseudo-terminal
                master_fd, slave_fd = pty.openpty()

                # 3. Launch VMD, attaching it to the PTY
                process = subprocess.Popen(
                    cmd,
                    stdin=slave_fd,
                    stdout=slave_fd,
                    stderr=slave_fd,
                    start_new_session=True,
                    close_fds=True,
                    cwd=vmd_cwd,
                    env=vmd_env
                )

                # 4. Close the slave end, VMD's process now owns it
                os.close(slave_fd)

                # 5. Start the background thread to read from our end
                reader_thread = threading.Thread(
                    target=pty_reader,
                    args=(master_fd,),
                    daemon=True  # Ensures thread exits when main script exits
                )
                reader_thread.start()

        except Exception as e:
            raise RuntimeError(f"VMD launch failed: {e}") from e

        return  # Prevent ChimeraX execution path

    if viewer_name != "chimerax":
        raise ValueError("viewer must be one of: 'chimerax', 'vmd', or 'molstar'.")

    #
    # =====================================================================
    # ------------------------- CHIMERAX BACKEND --------------------------
    # =====================================================================
    #

    chimera_commands = f"{open_commands} ; " + " ; ".join(default_commands)

    try:
        chx_exec = _which_chimerax()
    except FileNotFoundError as e:
        raise

    cmd = [chx_exec, '--cmd', chimera_commands]
    print(f"[flipbook] Launching: {cmd[0]}")

    try:
        # Use Popen to unblock the notebook, just like we do for VMD
        subprocess.Popen(cmd, start_new_session=True)
    except FileNotFoundError:
        raise FileNotFoundError("ChimeraX executable not found or not executable.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"ChimeraX exited with an error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred while running ChimeraX: {e}") from e


# -------------------------------------------------------------------------
# --------------------------- CLI INTERFACE -------------------------------
# -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Open PDB files in ChimeraX in numerical order and apply dynamic coloring and tiling.')
    parser.add_argument('directory', type=str, help='Path to the directory containing PDB files.')
    parser.add_argument('--palette', type=str, default='viridis', choices=COLOR_PALETTES.keys(),
                        help='Color palette to use for coloring the models.')
    parser.add_argument('--min_bfactor', type=float, default=None, help='Minimum B-factor value.')
    parser.add_argument('--max_bfactor', type=float, default=None, help='Maximum B-factor value.')
    parser.add_argument('--viewer', type=str, default='chimerax', choices=['chimerax', 'vmd', 'molstar'],
                        help='Choose visualization backend.')
    parser.add_argument('--molstar-output', type=str, default=None,
                        help='Output HTML path for the Molstar backend.')
    parser.add_argument('--molstar-manifest-output', type=str, default=None,
                        help='Output JSON manifest path for the Molstar backend.')
    parser.add_argument('--molstar-asset-mode', type=str, default='cdn', choices=['cdn', 'local', 'inline'],
                        help='How Molstar assets are loaded in generated HTML.')
    parser.add_argument('--molstar-height', type=int, default=720,
                        help='Notebook iframe height for the Molstar backend.')
    parser.add_argument('--molstar-camera-mode', type=str, default='orthographic',
                        choices=['orthographic', 'perspective'],
                        help='Camera projection mode for the Molstar backend.')
    parser.add_argument('--extra-commands', type=str, nargs='+', default=[],
                        help='Extra ChimeraX commands to run after the default commands.')
    args = parser.parse_args()

    try:
        run_flipbook(
            directory=args.directory,
            palette=args.palette,
            min_bfactor=args.min_bfactor,
            max_bfactor=args.max_bfactor,
            extra_commands=args.extra_commands,
            viewer=args.viewer,
            molstar_output=args.molstar_output,
            molstar_manifest_output=args.molstar_manifest_output,
            molstar_asset_mode=args.molstar_asset_mode,
            molstar_height=args.molstar_height,
            molstar_camera_mode=args.molstar_camera_mode
        )
    except Exception as e:
        print(f"[flipbook] ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
