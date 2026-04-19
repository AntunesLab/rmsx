import os
import shutil
from pathlib import Path


RMSX_OUTPUT_SENTINEL = ".rmsx_output_dir"
_RMSX_ALLOWED_TOP_LEVEL_NAMES = {"combined"}
_RMSX_BENIGN_TOP_LEVEL_FILES = {
    RMSX_OUTPUT_SENTINEL,
    ".DS_Store",
    "Thumbs.db",
}
_PREVIEW_LIMIT = 25


def _resolve_output_path(path_str):
    return Path(path_str).expanduser().resolve(strict=False)


def _path_is_within(parent_path, child_path):
    try:
        child_path.relative_to(parent_path)
        return True
    except ValueError:
        return False


def get_output_dir_safety_reason(output_dir, topology_file=None, trajectory_file=None):
    output_path = _resolve_output_path(output_dir)
    home_path = Path.home().resolve(strict=False)
    cwd_path = Path.cwd().resolve(strict=False)
    root_path = Path(output_path.anchor).resolve(strict=False) if output_path.anchor else None

    if root_path is not None and output_path == root_path:
        return (
            f"Refusing to overwrite '{output_path}' because it is a filesystem root. "
            "Choose a dedicated RMSX output directory."
        )

    if output_path == home_path:
        return (
            f"Refusing to overwrite '{output_path}' because it is your home directory. "
            "Choose a dedicated RMSX output directory."
        )

    if output_path == cwd_path:
        return (
            f"Refusing to overwrite '{output_path}' because it is the current working directory. "
            "Choose a dedicated RMSX output directory."
        )

    for input_file in (topology_file, trajectory_file):
        if not input_file:
            continue
        input_path = _resolve_output_path(input_file)
        input_parent = input_path.parent

        if output_path == input_parent:
            return (
                f"Refusing to overwrite '{output_path}' because it is the same directory as "
                f"the input file '{input_path.name}'. Choose a dedicated RMSX output directory."
            )

        if _path_is_within(output_path, input_path):
            return (
                f"Refusing to overwrite '{output_path}' because it contains the input file "
                f"'{input_path}'. Choose a dedicated RMSX output directory."
            )

    return None


def get_output_dir_overwrite_preview(output_dir, limit=_PREVIEW_LIMIT):
    output_path = Path(output_dir)

    if not output_path.exists() or not output_path.is_dir():
        return "No existing top-level entries would be deleted."

    entries = []
    for entry in sorted(output_path.iterdir(), key=lambda path: path.name.lower()):
        suffix = "/"
        if entry.is_file():
            suffix = ""
        elif entry.is_symlink():
            suffix = "@"
        entries.append(f"- {entry.name}{suffix}")

    if not entries:
        return "No existing top-level entries would be deleted."

    if len(entries) > limit:
        shown = entries[:limit]
        shown.append(f"- ... ({len(entries) - limit} more)")
        entries = shown

    return "Top-level entries that would be deleted:\n" + "\n".join(entries)


def is_rmsx_managed_output_dir(output_dir):
    output_path = Path(output_dir)

    if not output_path.exists():
        return True

    if not output_path.is_dir():
        return False

    for entry in output_path.iterdir():
        if entry.name in _RMSX_BENIGN_TOP_LEVEL_FILES:
            continue
        if entry.is_dir() and (
            entry.name in _RMSX_ALLOWED_TOP_LEVEL_NAMES
            or entry.name.startswith("chain_")
        ):
            continue
        return False

    return True


def write_output_dir_sentinel(output_dir):
    sentinel_path = Path(output_dir) / RMSX_OUTPUT_SENTINEL
    sentinel_path.write_text("managed_by=rmsx\n", encoding="utf-8")


def prepare_managed_output_dir(
    output_dir,
    overwrite=False,
    verbose=True,
    topology_file=None,
    trajectory_file=None,
):
    output_path = Path(output_dir)
    overwrite_preview = get_output_dir_overwrite_preview(output_path)
    safety_reason = get_output_dir_safety_reason(
        output_dir,
        topology_file=topology_file,
        trajectory_file=trajectory_file,
    )
    if safety_reason:
        if verbose and output_path.exists():
            print(overwrite_preview)
        raise RuntimeError(f"{safety_reason}\n{overwrite_preview}")

    if output_path.exists() and not output_path.is_dir():
        raise NotADirectoryError(
            f"Output path '{output_path}' exists but is not a directory."
        )

    if output_path.exists():
        if not is_rmsx_managed_output_dir(output_path):
            if verbose:
                print(overwrite_preview)
            raise RuntimeError(
                f"Refusing to overwrite '{output_path}' because it does not look like an "
                "RMSX-managed output directory. Allowed top-level contents are 'combined', "
                f"'chain_*', and the RMSX sentinel file.\n{overwrite_preview}"
            )

        if overwrite:
            if verbose:
                print(f"Clearing main output directory: {output_path}")
        else:
            response = input(
                f"The main directory '{output_path}' already exists. Overwrite? (y/n): "
            )
            if response.strip().lower() != "y":
                raise RuntimeError("User chose not to overwrite the main output directory.")
            if verbose:
                print(f"Clearing main output directory: {output_path}")

        for file in os.listdir(output_path):
            file_path = output_path / file
            try:
                if file_path.is_file() or file_path.is_symlink():
                    os.unlink(file_path)
                elif file_path.is_dir():
                    shutil.rmtree(file_path)
            except Exception as e:
                if verbose:
                    print(f"Failed to delete {file_path}. Reason: {e}")
    else:
        os.makedirs(output_path)
        if verbose:
            print(f"Created main output directory: {output_path}")

    write_output_dir_sentinel(output_path)
