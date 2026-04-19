import json
import os
import unittest
from pathlib import Path


class TestNotebookSafety(unittest.TestCase):
    def test_import_rmsx_without_core_deps(self) -> None:
        # Importing the top-level package should not eagerly import heavy deps
        # (e.g., pandas/MDAnalysis) so notebook users can at least access flipbook tools.
        import rmsx  # noqa: F401

    def test_run_flipbook_exposed_and_lightweight(self) -> None:
        from rmsx import run_flipbook

        self.assertEqual(run_flipbook.__module__, "rmsx.flipbook")

    def test_run_flipbook_raises_not_sysexit(self) -> None:
        from rmsx import run_flipbook

        with self.assertRaises(NotADirectoryError):
            run_flipbook("/this/does/not/exist", viewer="chimerax")

    def test_vmd_scripts_packaged_and_loader_is_cwd_independent(self) -> None:
        import rmsx.vmd_scripts as vmd_scripts

        scripts_dir = Path(vmd_scripts.__file__).resolve().parent
        wait_to_load = scripts_dir / "wait_to_load.tcl"
        grid_script = scripts_dir / "grid_color_scale_centered_xaxis_hotkeys.tcl"

        self.assertTrue(wait_to_load.is_file(), "wait_to_load.tcl should be present next to vmd_scripts package")
        self.assertTrue(grid_script.is_file(), "grid TCL script should be present next to vmd_scripts package")

        text = wait_to_load.read_text(encoding="utf-8")
        self.assertIn("[info script]", text)
        self.assertIn("after idle", text)
        self.assertIn("VMDMODULATENEWTUBE", text)

    def test_flipbook_can_find_loader_script_on_disk(self) -> None:
        # When installed normally (wheel unpacked), this path should exist.
        # In this repo checkout, it should also exist.
        import rmsx.flipbook as flipbook

        loader = os.path.join(os.path.dirname(flipbook.__file__), "vmd_scripts", "wait_to_load.tcl")
        self.assertTrue(os.path.exists(loader))

    def test_packaged_demo_inputs_exist_next_to_module(self) -> None:
        import rmsx

        pkg_dir = Path(rmsx.__file__).resolve().parent
        demo_dir = pkg_dir / "test_files"

        self.assertTrue((demo_dir / "1UBQ.pdb").is_file())
        self.assertTrue((demo_dir / "mon_sys.dcd").is_file())
        self.assertTrue((demo_dir / "protease_backbone.pdb").is_file())
        self.assertTrue((demo_dir / "short_protease_backbone.dcd").is_file())

    def test_notebook_uses_packaged_demo_inputs_and_external_output_root(self) -> None:
        notebook_path = Path("/Users/finn/Documents/GitHub/rmsx/RMSX_FlipBook_Quickstart.ipynb")
        nb = json.loads(notebook_path.read_text(encoding="utf-8"))
        source = "\n".join("".join(cell.get("source", [])) for cell in nb["cells"])

        self.assertIn('pkg_dir / "test_files"', source)
        self.assertIn('demo_output_root = Path.cwd() / "rmsx_demo_outputs"', source)
        self.assertIn('output_dir = (demo_output_root / "example_uqb").as_posix()', source)
        self.assertIn('output_dir_multi = (demo_output_root / "protease").as_posix()', source)
        self.assertNotIn('output_dir = (test_dir / "example_uqb").as_posix()', source)
        self.assertNotIn('output_dir_multi = (test_dir / "protease").as_posix()', source)

    def test_notebook_includes_advanced_masking_example(self) -> None:
        notebook_path = Path("/Users/finn/Documents/GitHub/rmsx/RMSX_FlipBook_Quickstart.ipynb")
        nb = json.loads(notebook_path.read_text(encoding="utf-8"))
        source = "\n".join("".join(cell.get("source", [])) for cell in nb["cells"])

        self.assertIn("## 7) Advanced Feature: Masking Selected Residues", source)
        self.assertIn("protease_active_site_mask = [", source)
        self.assertIn("segid A and resid 45:55", source)
        self.assertIn("segid B and resid 45:55", source)
        self.assertIn('output_dir_mask_demo = Path.cwd() / "rmsx_demo_outputs" / "protease_mask_example"', source)
        self.assertIn("all_chain_rmsx(", source)
        self.assertIn("mask=protease_active_site_mask", source)
        self.assertIn("resave the topology/trajectory without that chain or region using **MDAnalysis**", source)


if __name__ == "__main__":
    unittest.main()
