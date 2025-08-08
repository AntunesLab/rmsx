#!/usr/bin/env python3
import argparse
from rmsx import run_rmsx

def main():
    p = argparse.ArgumentParser(
        description="RMSX Trajectory Analysis Tool (MDAnalysis + R plotting + ChimeraX)"
    )
    p.add_argument("--topology", required=True, help="Topology file (PDB/GRO/PSF/etc.)")
    p.add_argument("--trajectory", required=True, help="Trajectory file (DCD/XTC/TRR/etc.)")
    p.add_argument("--output-dir", default=None, help="Output directory (default: <toponame>_rmsx)")
    # choose ONE: --num-slices OR --slice-size
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--num-slices", type=int, help="Number of equal time slices")
    g.add_argument("--slice-size", type=int, help="Slice size in frames")
    p.add_argument("--chain", dest="chain_sele", default=None, help="Chain/segid to analyze (e.g., A or 7)")
    p.add_argument("--palette", default="viridis", help="R/viridis palette (e.g., viridis, plasma, mako, turbo)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite an existing output directory")
    p.add_argument("--start-frame", type=int, default=0, help="First frame index to analyze (default: 0)")
    p.add_argument("--end-frame", type=int, default=None, help="Last frame index to analyze (inclusive)")
    p.add_argument("--no-plot", action="store_true", help="Skip calling the R plotting script")
    p.add_argument("--triple", action="store_true", help="Include RMSD & RMSF plots next to the heatmap")
    p.add_argument("--interpolate", action="store_true", help="Interpolate heatmap between time slices")
    p.add_argument("--log-transform", action="store_true", help="Log-transform RMSX before writing PDB B-factors")
    p.add_argument("--rscript", default="Rscript", help="Path to Rscript (default: Rscript on PATH)")
    p.add_argument("--verbose", action="store_true", help="Verbose logs")

    args = p.parse_args()

    run_rmsx(
        topology_file=args.topology,
        trajectory_file=args.trajectory,
        output_dir=args.output_dir,
        num_slices=args.num_slices,
        slice_size=args.slice_size,
        rscript_executable=args.rscript,
        verbose=args.verbose,
        interpolate=args.interpolate,
        triple=args.triple,
        chain_sele=args.chain_sele,
        overwrite=args.overwrite,
        palette=args.palette,
        start_frame=args.start_frame,
        end_frame=args.end_frame,
        make_plot=(not args.no_plot),
        log_transform=args.log_transform
    )

if __name__ == "__main__":
    main()
