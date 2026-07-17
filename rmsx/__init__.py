from __future__ import annotations

from importlib import import_module
from typing import Any

__all__ = [
    "run_rmsx",
    "combine_pdb_files",
    "all_chain_rmsx",
    "run_rmsx_flipbook",
    "run_flipbook",
    "run_shift_map",
    "all_chain_shift_map",
    "run_shift_flipbook",
    "MolstarFlipbookResult",
    "build_molstar_manifest",
    "write_molstar_flipbook",
]

_CORE_EXPORTS = {
    "run_rmsx",
    "combine_pdb_files",
    "all_chain_rmsx",
    "run_rmsx_flipbook",
    "run_shift_map",
    "all_chain_shift_map",
    "run_shift_flipbook",
}


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        core = import_module(".core", __name__)
        return getattr(core, name)
    if name == "run_flipbook":
        flipbook = import_module(".flipbook", __name__)
        return getattr(flipbook, name)
    if name in {"MolstarFlipbookResult", "build_molstar_manifest", "write_molstar_flipbook"}:
        molstar_viewer = import_module(".molstar_viewer", __name__)
        return getattr(molstar_viewer, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    # Avoid importing heavy dependencies during introspection/test discovery.
    return sorted(set(globals().keys()))
