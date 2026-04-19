import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest import mock

import rmsx.flipbook as flipbook
from rmsx.vmd_scripts.vmd_finder import find_vmd_executable


def _write_minimal_pdb(path: Path, chain_id: str, residue_ids) -> None:
    lines = []
    for serial_number, residue_id in enumerate(residue_ids, start=1):
        lines.append(
            f"ATOM  {serial_number:5d}  CA  ALA {chain_id}{residue_id:4d}    "
            f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}{10.00:6.2f}           C\n"
        )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


class TestFlipbookMasking(unittest.TestCase):
    def test_build_chimerax_mask_commands_maps_numeric_chain_to_x(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            pdb_path = tmp_path / "slice_1_first_frame.pdb"
            _write_minimal_pdb(pdb_path, "7", [10, 11, 12])
            (tmp_path / flipbook.MASK_METADATA_FILENAME).write_text(
                "ResidueID,ChainID,Masked\n10,7,True\n11,7,True\n12,7,True\n",
                encoding="utf-8",
            )

            commands = flipbook._build_chimerax_mask_commands(tmp_path, [str(pdb_path)])

            self.assertEqual(commands, ["transparency /X:10-12 70 target r"])

    def test_run_flipbook_appends_mask_commands_for_chimerax(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            pdb_path = tmp_path / "slice_1_first_frame.pdb"
            _write_minimal_pdb(pdb_path, "A", [1, 2, 3])
            (tmp_path / flipbook.MASK_METADATA_FILENAME).write_text(
                "ResidueID,ChainID,Masked\n1,A,True\n2,A,False\n3,A,True\n",
                encoding="utf-8",
            )

            with mock.patch.object(flipbook, "_which_chimerax", return_value="/usr/bin/chimerax"), \
                 mock.patch.object(flipbook.subprocess, "Popen") as popen_mock:
                flipbook.run_flipbook(str(tmp_path), viewer="chimerax")

            call_args, _ = popen_mock.call_args
            called_command = call_args[0]
            command_text = called_command[2]
            self.assertIn("transparency /A:1,3 70 target r", command_text)

    def test_run_flipbook_sets_vmd_mask_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            pdb_path = tmp_path / "slice_1_first_frame.pdb"
            _write_minimal_pdb(pdb_path, "A", [1, 2])
            (tmp_path / flipbook.MASK_METADATA_FILENAME).write_text(
                "ResidueID,ChainID,Masked\n1,A,True\n2,A,False\n",
                encoding="utf-8",
            )

            with mock.patch.object(flipbook, "find_vmd_executable", return_value="/usr/bin/vmd"), \
                 mock.patch.object(flipbook, "pty", None), \
                 mock.patch.object(flipbook.subprocess, "Popen") as popen_mock:
                flipbook.run_flipbook(str(tmp_path), viewer="vmd")

            _, call_kwargs = popen_mock.call_args
            environment = call_kwargs["env"]
            self.assertEqual(environment["RMSX_VMD_MASK_FILE"], str(tmp_path / flipbook.MASK_METADATA_FILENAME))
            self.assertEqual(environment["RMSX_VMD_MASK_OPACITY"], str(flipbook.VMD_MASK_OPACITY))

    def test_vmd_mask_tcl_creates_transparent_rep_when_vmd_available(self) -> None:
        try:
            vmd_exec = find_vmd_executable()
        except FileNotFoundError:
            self.skipTest("VMD executable not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            pdb1 = tmp_path / "slice_1_first_frame.pdb"
            pdb2 = tmp_path / "slice_2_first_frame.pdb"
            _write_minimal_pdb(pdb1, "A", [1, 2, 3])
            _write_minimal_pdb(pdb2, "A", [1, 2, 3])
            mask_path = tmp_path / flipbook.MASK_METADATA_FILENAME
            mask_path.write_text(
                "ResidueID,ChainID,Masked\n1,A,True\n2,A,False\n3,A,True\n",
                encoding="utf-8",
            )

            wrapper = tmp_path / "vmd_mask_smoke.tcl"
            wrapper.write_text(textwrap.dedent("""
                source /Users/finn/Documents/GitHub/rmsx/rmsx/vmd_scripts/wait_to_load.tcl
                run_main_script_after_delay
                set mid [lindex $::molList 0]
                puts "SMOKE numreps=[molinfo $mid get numreps]"
                puts "SMOKE rep0_selection=[molinfo $mid get "{selection 0}"]"
                puts "SMOKE rep1_selection=[molinfo $mid get "{selection 1}"]"
                puts "SMOKE rep1_material=[molinfo $mid get "{material 1}"]"
                quit
            """), encoding="utf-8")

            env = os.environ.copy()
            env["RMSX_VMD_MAIN"] = "/Users/finn/Documents/GitHub/rmsx/rmsx/vmd_scripts/grid_color_scale_centered_xaxis_hotkeys.tcl"
            env["RMSX_VMD_MASK_FILE"] = str(mask_path)
            env["RMSX_VMD_MASK_OPACITY"] = str(flipbook.VMD_MASK_OPACITY)
            result = subprocess.run(
                [
                    vmd_exec,
                    "-dispdev", "text",
                    "-e", str(wrapper),
                    "-args",
                    str(pdb1),
                    str(pdb2),
                    "viridis",
                ],
                capture_output=True,
                text=True,
                env=env,
                check=False,
            )

            output = result.stdout + result.stderr
            self.assertIn("Masked transparency applied:", output)
            self.assertIn("SMOKE numreps=2", output)
            self.assertIn("SMOKE rep1_material=Transparent", output)
            self.assertIn("resid 1 3", output)


if __name__ == "__main__":
    unittest.main()
