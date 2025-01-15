# rmsx/rmsx/__init__.py

from .rmsx import run_rmsx, combine_pdb_files, plot_rmsx_surface, all_chain_rmsx, run_rmsx_flipbook
from .flipbook import run_flipbook  # Add this line to import run_flipbook


# updated to include flipbook option
# Update __all__ to include run_flipbook
__all__ = [
    'run_rmsx',
    'combine_pdb_files',
    'plot_rmsx_surface',
    'all_chain_rmsx',
    'run_rmsx_flipbook',
    'run_flipbook'
]
