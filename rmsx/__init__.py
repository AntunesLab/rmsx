# rmsx/rmsx/__init__.py

from .rmsx import run_rmsx, combine_pdb_files, all_chain_rmsx, run_rmsx_flipbook, run_shift_map, all_chain_shift_map, run_shift_flipbook
from .flipbook import run_flipbook  # Add this line to import run_flipbook


# updated to include flipbook option
# Update __all__ to include run_flipbook
__all__ = [
    'run_rmsx',
    'combine_pdb_files', #     removed 'plot_rmsx_surface' at some point
    'all_chain_rmsx',
    'run_rmsx_flipbook',
    'run_flipbook',
    'run_shift_map',
    'all_chain_shift_map',
    'run_shift_flipbook'

]
