#!/usr/bin/env python3
import argparse
from importlib.metadata import version, PackageNotFoundError
from .core import run_rmsx  # import directly from the module to avoid any circulars

def _pkg_version() -> str:
    try:
        return version("rmsx")
    except PackageNotFoundError:
        # fallback if running from source without an installed dist-info
        try:
            from . import __version__
            return __version__  # type: ignore[attr-defined]
        except Exception:
            return "0.0.0+local"

def build_parser(prog: str = "rmsx") -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=prog,
        description="RMSX Trajectory Analysis Tool",
    )
    p.add_argument("psf_file", help="Path to the PSF file")
    p.add_argument("dcd_file", help="Path to the DCD file")
    p.add_argument("pdb_file", help="Path to the PDB file")
    p.add_argument("--output_dir", default=None, help="Output directory")
    # NOTE: use float if run_rmsx expects fractional slice sizes; keep int if not.
    p.add_argument("--slice_size", type=int, default=5, help="Slice size for trajectory processing")
    p.add_argument("--verbose", action="store_true", help="Enable verbose output")
    p.add_argument("--interpolate", action="store_true", help="Enable interpolation in plots")
    p.add_argument("--triple", action="store_true", help="Enable triple plotting")
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s { _pkg_version() }")
    return p

def main(argv=None, prog: str = "rmsx") -> None:
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)
    run_rmsx(
        psf_file=args.psf_file,
        dcd_file=args.dcd_file,
        pdb_file=args.pdb_file,
        output_dir=args.output_dir,
        slice_size=args.slice_size,
        verbose=args.verbose,
        interpolate=args.interpolate,
        triple=args.triple,
    )

if __name__ == "__main__":
    main(prog="rmsx")

