import tempfile
import unittest
from pathlib import Path

from rmsx.output_safety import (
    RMSX_OUTPUT_SENTINEL,
    get_output_dir_safety_reason,
    is_rmsx_managed_output_dir,
    prepare_managed_output_dir,
)


class TestOutputDirSafety(unittest.TestCase):
    def test_flags_input_directory_as_unsafe(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            input_dir = base / "inputs"
            input_dir.mkdir()
            topology = input_dir / "prot.pdb"
            trajectory = input_dir / "traj.dcd"
            topology.write_text("", encoding="utf-8")
            trajectory.write_text("", encoding="utf-8")

            reason = get_output_dir_safety_reason(
                input_dir,
                topology_file=topology,
                trajectory_file=trajectory,
            )

            self.assertIsNotNone(reason)
            self.assertIn("same directory as the input file", reason)

    def test_flags_parent_of_input_directory_as_unsafe(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            input_dir = base / "inputs"
            input_dir.mkdir()
            topology = input_dir / "prot.pdb"
            trajectory = input_dir / "traj.dcd"
            topology.write_text("", encoding="utf-8")
            trajectory.write_text("", encoding="utf-8")

            reason = get_output_dir_safety_reason(
                base,
                topology_file=topology,
                trajectory_file=trajectory,
            )

            self.assertIsNotNone(reason)
            self.assertIn("contains the input file", reason)

    def test_recognizes_managed_directory_layout(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "chain_A_rmsx").mkdir()
            (output_dir / "chain_B_lddtmap").mkdir()
            (output_dir / "combined").mkdir()
            (output_dir / ".DS_Store").write_text("", encoding="utf-8")

            self.assertTrue(is_rmsx_managed_output_dir(output_dir))

    def test_rejects_unmanaged_directory_layout(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "chain_A_rmsx").mkdir()
            (output_dir / "notes").mkdir()

            self.assertFalse(is_rmsx_managed_output_dir(output_dir))

    def test_prepare_managed_output_dir_refuses_unmanaged_existing_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "notes").mkdir()
            (output_dir / "results.txt").write_text("", encoding="utf-8")

            with self.assertRaises(RuntimeError) as ctx:
                prepare_managed_output_dir(output_dir, overwrite=True, verbose=False)

            message = str(ctx.exception)
            self.assertIn("does not look like an RMSX-managed output directory", message)
            self.assertIn("Top-level entries that would be deleted:", message)
            self.assertIn("- notes/", message)
            self.assertIn("- results.txt", message)

    def test_prepare_managed_output_dir_previews_dangerous_input_directory_contents(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir = Path(tmpdir) / "inputs"
            input_dir.mkdir()
            topology = input_dir / "prot.pdb"
            trajectory = input_dir / "traj.dcd"
            extra_file = input_dir / "notes.txt"
            topology.write_text("", encoding="utf-8")
            trajectory.write_text("", encoding="utf-8")
            extra_file.write_text("", encoding="utf-8")

            with self.assertRaises(RuntimeError) as ctx:
                prepare_managed_output_dir(
                    input_dir,
                    overwrite=True,
                    verbose=False,
                    topology_file=topology,
                    trajectory_file=trajectory,
                )

            message = str(ctx.exception)
            self.assertIn("same directory as the input file", message)
            self.assertIn("Top-level entries that would be deleted:", message)
            self.assertIn("- prot.pdb", message)
            self.assertIn("- traj.dcd", message)
            self.assertIn("- notes.txt", message)

    def test_prepare_managed_output_dir_clears_managed_folder_and_rewrites_sentinel(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            input_dir = base / "inputs"
            output_dir = base / "results"
            input_dir.mkdir()
            output_dir.mkdir()

            topology = input_dir / "prot.pdb"
            trajectory = input_dir / "traj.dcd"
            topology.write_text("", encoding="utf-8")
            trajectory.write_text("", encoding="utf-8")

            (output_dir / "chain_A_rmsx").mkdir()
            (output_dir / "combined").mkdir()

            prepare_managed_output_dir(
                output_dir,
                overwrite=True,
                verbose=False,
                topology_file=topology,
                trajectory_file=trajectory,
            )

            self.assertTrue((output_dir / RMSX_OUTPUT_SENTINEL).is_file())
            remaining = {path.name for path in output_dir.iterdir()}
            self.assertEqual(remaining, {RMSX_OUTPUT_SENTINEL})


if __name__ == "__main__":
    unittest.main()
