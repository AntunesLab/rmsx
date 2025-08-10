#!/usr/bin/env python3
import argparse, sys
from importlib.metadata import version, PackageNotFoundError
from .core import run_rmsx

def _pkg_version() -> str:
    try:
        return version("rmsx")
    except PackageNotFoundError:
        try:
            from . import __version__  # optional fallback
            return __version__         # type: ignore[attr-defined]
        except Exception:
            return "0.0.0+local"

def build_parser(prog: str = "rmsx") -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog=prog, description="RMSX Trajectory Analysis Tool")

    # required (match run_rmsx signature: topology_file, trajectory_file)
    p.add_argument("psf_file", help="Topology file (e.g., PSF/PDB/PRMTOP)")
    p.add_argument("dcd_file", help="Trajectory file (e.g., DCD/XTC/TRR/…)")

    # optional legacy 3rd positional (ignored)
    p.add_argument("pdb_file", nargs="?", help=argparse.SUPPRESS)

    # common options (None means 'only pass if user set it')
    p.add_argument("--output_dir", default=None, help="Output directory")
    p.add_argument("--slice_size", type=int, default=None, help="Slice size (frames)")
    p.add_argument("--num_slices", type=int, default=None, help="Number of slices (overrides slice_size)")
    p.add_argument("--chain", dest="chain_sele", default=None, help="Chain selection (e.g., 'A' or 'A,B')")
    p.add_argument("--palette", default=None, help="Palette name (default: viridis)")
    p.add_argument("--start_frame", type=int, default=None, help="Start frame index")
    p.add_argument("--end_frame", type=int, default=None, help="End frame index")
    p.add_argument("--analysis_type", choices=["protein","dna","rna","generic"], default=None)
    p.add_argument("--summary_n", type=int, default=None)
    p.add_argument("--manual_length_ns", type=float, default=None)
    p.add_argument("--custom_fill_label", default=None)
    p.add_argument("--rscript", dest="rscript_executable", default=None, help="Path to Rscript")

    # booleans as tri-state (None unless explicitly set)
    vb = p.add_mutually_exclusive_group()
    vb.add_argument("--verbose", dest="verbose", action="store_true")
    vb.add_argument("--quiet",   dest="verbose", action="store_false")

    ib = p.add_mutually_exclusive_group()
    ib.add_argument("--interpolate",    dest="interpolate", action="store_true")
    ib.add_argument("--no-interpolate", dest="interpolate", action="store_false")

    p.add_argument("--triple", action="store_true", default=None)
    p.add_argument("--log_transform", action="store_true", default=None)
    p.add_argument("--no-plot", dest="make_plot", action="store_false", default=None)
    p.add_argument("--overwrite", action="store_true", default=None)

    # make unspecified flags stay as None
    p.set_defaults(verbose=None, interpolate=None)

    p.add_argument("-V", "--version", action="version", version=f"%(prog)s { _pkg_version() }")
    return p

def main(argv=None, prog: str = "rmsx") -> None:
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    # Gentle note if a legacy 3rd positional was provided
    if getattr(args, "pdb_file", None):
        print("[rmsx] Note: ignoring third positional (pdb_file); API uses topology & trajectory only.",
              file=sys.stderr)

    # Only pass kwargs the user actually set (avoid clobbering core defaults)
    kw = {}
    if args.output_dir is not None:        kw["output_dir"] = args.output_dir
    if args.num_slices is not None:        kw["num_slices"] = args.num_slices
    if args.slice_size is not None:        kw["slice_size"] = args.slice_size
    if args.rscript_executable is not None:kw["rscript_executable"] = args.rscript_executable
    if args.verbose is not None:           kw["verbose"] = args.verbose
    if args.interpolate is not None:       kw["interpolate"] = args.interpolate
    if args.triple is not None:            kw["triple"] = args.triple
    if args.chain_sele is not None:        kw["chain_sele"] = args.chain_sele
    if args.overwrite is not None:         kw["overwrite"] = args.overwrite
    if args.palette is not None:           kw["palette"] = args.palette
    if args.start_frame is not None:       kw["start_frame"] = args.start_frame
    if args.end_frame is not None:         kw["end_frame"] = args.end_frame
    if args.make_plot is not None:         kw["make_plot"] = args.make_plot
    if args.analysis_type is not None:     kw["analysis_type"] = args.analysis_type
    if args.summary_n is not None:         kw["summary_n"] = args.summary_n
    if args.manual_length_ns is not None:  kw["manual_length_ns"] = args.manual_length_ns
    if args.log_transform is not None:     kw["log_transform"] = args.log_transform
    if args.custom_fill_label is not None: kw["custom_fill_label"] = args.custom_fill_label

    # run_rmsx(topology_file, trajectory_file, **options)
    run_rmsx(args.psf_file, args.dcd_file, **kw)

if __name__ == "__main__":
    main(prog="rmsx")

