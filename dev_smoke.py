#!/usr/bin/env python3
from pathlib import Path
import inspect
import traceback
import sys

import rmsx

ROOT = Path(__file__).resolve().parent
TF   = ROOT / "test_files"
OUT  = TF / "_smoke_out"
OUT.mkdir(exist_ok=True)

# Inputs
PSF = TF / "1UBQ.psf"
DCD = TF / "mon_sys.dcd"
PDB = TF / "1UBQ.pdb"

# A small directory with lots of *first_frame.pdb files
EXAMPLE_CHAIN_DIR = TF / "example1" / "chain_7_rmsx"
EXAMPLE_CHAIN_DIR.mkdir(parents=True, exist_ok=True)  # it should already exist

def build_kwargs(func, extra=None):
    """Return kwargs subset that func accepts, avoiding clobbering defaults with None."""
    extra = extra or {}
    base = dict(
        output_dir=OUT / func.__name__,
        slice_size=5,
        num_slices=None,
        verbose=True,
        interpolate=False,
        triple=False,
        overwrite=True,              # filtered out if func doesn't accept it
        start_frame=0,
        end_frame=40,
        make_plot=False,
        analysis_type="protein",
        summary_n=1,
        palette="viridis",
        manual_length_ns=None,
        log_transform=False,
        custom_fill_label=None,
        chain_sele=None,
        rscript_executable=None,
        launch=False,                # avoid launching ChimeraX if supported
        # convenience alternates some funcs might expose:
        topology=PSF, trajectory=DCD, structure=PDB,
        pdb_file=PDB,
        directory=str(EXAMPLE_CHAIN_DIR),
        chain_dirs=[str(EXAMPLE_CHAIN_DIR)],
        output_file=str((OUT / func.__name__) / "combined.pdb"),
    )
    base.update(extra or {})
    sig = inspect.signature(func)
    names = set(sig.parameters.keys())
    return {k: v for k, v in base.items() if k in names and v is not None}

def call(func):
    sig = inspect.signature(func)
    params = list(sig.parameters.keys())

    kwargs = build_kwargs(func)

    # If the first two params look like topology/trajectory, feed them POSITIONALLY
    # and make sure we don't also pass them as kwargs.
    first_two = params[:2]
    if first_two == ["topology_file", "trajectory_file"] or first_two == ["topology", "trajectory"]:
        for k in ("topology_file", "trajectory_file", "topology", "trajectory"):
            kwargs.pop(k, None)
        return func(PSF, DCD, **kwargs)

    # run_flipbook likely needs a directory; keep headless if possible
    if func.__name__ in {"run_flipbook"}:
        if "directory" in params and "directory" not in kwargs:
            kwargs["directory"] = str(EXAMPLE_CHAIN_DIR)
        if "launch" in params and "launch" not in kwargs:
            kwargs["launch"] = False

    # combine_pdb_files: support both new and old signatures
    if func.__name__ == "combine_pdb_files":
        # ensure outputs/inputs exist
        out_dir = OUT / func.__name__
        out_dir.mkdir(parents=True, exist_ok=True)
        chain_dirs_val = [str(EXAMPLE_CHAIN_DIR)]
        combined_dir_val = str(out_dir / "combined_dir")
        Path(combined_dir_val).mkdir(parents=True, exist_ok=True)

        if "chain_dirs" in params and "combined_dir" in params:
            return func(chain_dirs=chain_dirs_val, combined_dir=combined_dir_val)
        elif "chain_dirs" in params and "output_file" in params:
            return func(chain_dirs=chain_dirs_val, output_file=str(out_dir / "combined.pdb"))
        elif len(params) >= 2:
            # positional fallback: (chain_dirs, combined_dir)
            return func(chain_dirs_val, combined_dir_val)
        else:
            # last resort: try kwargs we already built (older bespoke signatures)
            return func(**kwargs)

    # Flipbook variants: prefer not launching during smoke
    if func.__name__ in {"run_rmsx_flipbook", "run_shift_flipbook"}:
        if "launch" in params and "launch" not in kwargs:
            kwargs["launch"] = False

    return func(**kwargs)

def run(name):
    f = getattr(rmsx, name, None)
    if f is None:
        print(f"SKIP {name}: not exported")
        return True
    try:
        print(f"→ {name} …", end=" ", flush=True)
        call(f)
        print("PASS")
        return True
    except Exception:
        print("FAIL")
        traceback.print_exc()
        return False

def main():
    print("Using:", sys.executable)
    print("rmsx from:", rmsx.__file__)
    print("Inputs:", PSF.exists(), DCD.exists(), PDB.exists(), EXAMPLE_CHAIN_DIR.exists())
    passed = True

    # core + variants
    for fn in [
        "run_rmsx",
        "all_chain_rmsx",
        "run_shift_map",
        "all_chain_shift_map",
        "run_rmsx_flipbook",
        "run_shift_flipbook",
        "run_flipbook",
        "combine_pdb_files",
    ]:
        passed &= run(fn)

    print("\nRESULT:", "ALL PASS" if passed else "FAILURES PRESENT")
    print("Output dir:", OUT)

if __name__ == "__main__":
    main()

