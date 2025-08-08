#!python

import argparse
from rmsx.rmsx import run_rmsx

def main():
    parser = argparse.ArgumentParser(description='RMSX Trajectory Analysis Tool')
    parser.add_argument('psf_file', help='Path to the PSF file')
    parser.add_argument('dcd_file', help='Path to the DCD file')
    parser.add_argument('pdb_file', help='Path to the PDB file')
    parser.add_argument('--output_dir', help='Output directory', default=None)
    parser.add_argument('--slice_size', type=int, default=5, help='Slice size for trajectory processing')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--interpolate', action='store_true', help='Enable interpolation in plots')
    parser.add_argument('--triple', action='store_true', help='Enable triple plotting')
    args = parser.parse_args()

    run_rmsx(
        psf_file=args.psf_file,
        dcd_file=args.dcd_file,
        pdb_file=args.pdb_file,
        output_dir=args.output_dir,
        slice_size=args.slice_size,
        verbose=args.verbose,
        interpolate=args.interpolate,
        triple=args.triple
    )

if __name__ == '__main__':
    main()

