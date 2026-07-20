import json
import tempfile
import unittest
from pathlib import Path

import rmsx.flipbook as flipbook
from rmsx.molstar_viewer import ASSET_DIR, build_molstar_manifest


def _atom_line(serial, atom_name, chain_id, residue_id, bfactor, x, y=0.0, z=0.0, segid=""):
    return (
        f"ATOM  {serial:5d} {atom_name:^4s} ALA {chain_id:1s}{residue_id:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}"
        f"      {segid:<4s}{atom_name[0]:>2s}  \n"
    )


def _write_slice(path: Path, chain_id: str, bfactors, segid="") -> None:
    lines = []
    serial = 1
    for residue_id, bfactor in enumerate(bfactors, start=1):
        lines.append(_atom_line(serial, "N", chain_id, residue_id, bfactor, float(residue_id), segid=segid))
        serial += 1
        lines.append(_atom_line(serial, "CA", chain_id, residue_id, bfactor, float(residue_id) + 0.25, segid=segid))
        serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


class TestMolstarViewer(unittest.TestCase):
    def test_slices_are_rigidly_aligned_to_first_slice(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            reference = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
            mobile = [(10.0 - y, 5.0 + x, 3.0 + z) for x, y, z in reference]
            for slice_index, coordinates in enumerate((reference, mobile), start=1):
                lines = [
                    _atom_line(index, "CA", "A", index, float(slice_index), x, y, z)
                    for index, (x, y, z) in enumerate(coordinates, start=1)
                ]
                source_text = "".join(lines) + "END\n"
                (tmp_path / f"slice_{slice_index}_first_frame.pdb").write_text(source_text, encoding="utf-8")

            manifest = build_molstar_manifest(tmp_path)
            aligned_lines = [
                line for line in manifest["slices"][1]["pdb"].splitlines() if line.startswith("ATOM")
            ]
            aligned = [
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                for line in aligned_lines
            ]

            self.assertEqual(manifest["alignmentModel"]["alignedSlices"], 2)
            self.assertEqual(manifest["alignmentModel"]["skippedSlices"], [])
            self.assertEqual(
                (tmp_path / "slice_2_first_frame.pdb").read_text(encoding="utf-8"),
                "".join(
                    _atom_line(index, "CA", "A", index, 2.0, x, y, z)
                    for index, (x, y, z) in enumerate(mobile, start=1)
                )
                + "END\n",
            )
            for observed, expected in zip(aligned, reference):
                for observed_value, expected_value in zip(observed, expected):
                    self.assertAlmostEqual(observed_value, expected_value, places=3)

    def test_default_rotation_uses_a_portrait_principal_axis(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            coordinates = [(0.0, 0.0, 0.0), (2.0, 2.0, 0.2), (4.0, 4.0, -0.1), (6.0, 6.0, 0.1)]
            lines = [
                _atom_line(index, "CA", "A", index, 1.0, x, y, z)
                for index, (x, y, z) in enumerate(coordinates, start=1)
            ]
            (tmp_path / "slice_1_first_frame.pdb").write_text("".join(lines) + "END\n", encoding="utf-8")

            manifest = build_molstar_manifest(tmp_path)
            rotation = manifest["rotationModel"]
            matrix = rotation["defaultRotationMatrix"]
            diagonal = (1.0, 1.0, 0.0)
            transformed = [sum(matrix[row][column] * diagonal[column] for column in range(3)) for row in range(3)]

            self.assertEqual(rotation["defaultRotationSource"], "principal-axis portrait orientation")
            self.assertEqual(rotation["orientation"], "portrait")
            self.assertEqual(rotation["orientationSelection"], "alpha carbons across aligned slices")
            self.assertAlmostEqual(transformed[0], 0.0, places=2)
            self.assertGreater(transformed[1], 1.4)
            self.assertAlmostEqual(transformed[2], 0.0, places=2)

    def test_run_flipbook_molstar_writes_manifest_and_notebook_html(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            _write_slice(tmp_path / "slice_1_first_frame.pdb", "A", [0.25, 0.50])
            _write_slice(tmp_path / "slice_2_first_frame.pdb", "A", [1.25, 2.50])

            result = flipbook.run_flipbook(
                str(tmp_path),
                viewer="molstar",
                palette="turbo",
                molstar_height=640,
            )

            self.assertTrue(result.html_path.is_file())
            self.assertTrue(result.manifest_path.is_file())
            self.assertIn('title="RMSX Molstar Flipbook"', result._repr_html_())
            self.assertIn('height="640"', result._repr_html_())
            self.assertIn("cdn.jsdelivr.net/npm/molstar@5.4.2", result.html)
            self.assertIn("dataset.incoming", result.html)

            manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["schemaVersion"], "flipbook-molstar-viewer/v1")
            self.assertEqual(manifest["palette"]["name"], "turbo")
            self.assertEqual(len(manifest["slices"]), 2)
            self.assertEqual(manifest["domain"], {"min": 0.25, "max": 2.5})
            self.assertEqual(manifest["summaries"]["slice_2.dcd"]["maxResidue"], "2")
            self.assertEqual(manifest["residues"][1]["values"]["slice_2.dcd"], 2.5)
            self.assertEqual(manifest["molstarRenderStyle"]["cameraMode"], "orthographic")
            self.assertEqual(manifest["flipbookReference"]["minimumSpacingFactor"], 0.0)
            self.assertEqual(manifest["flipbookReference"]["maximumSpacingFactor"], 1.5)
            self.assertEqual(manifest["flipbookReference"]["defaultSpacingMode"], "auto")
            self.assertIsNone(manifest["flipbookReference"]["defaultSpacingFactor"])
            self.assertEqual(manifest["flipbookReference"]["tilePaddingFactor"], 1.0)
            self.assertEqual(manifest["keyboardShortcuts"]["spacingStep"], 0.01)

            perspective_result = flipbook.run_flipbook(
                str(tmp_path),
                viewer="molstar",
                palette="turbo",
                molstar_output=tmp_path / "perspective.html",
                molstar_manifest_output=tmp_path / "perspective.json",
                molstar_camera_mode="perspective",
            )
            self.assertEqual(
                perspective_result.manifest["molstarRenderStyle"]["cameraMode"],
                "perspective",
            )

            with self.assertRaisesRegex(ValueError, "camera_mode"):
                build_molstar_manifest(tmp_path, camera_mode="fisheye")
            with self.assertRaisesRegex(ValueError, "orientation"):
                build_molstar_manifest(tmp_path, orientation="diagonal")

            fixed_spacing = build_molstar_manifest(
                tmp_path,
                spacing_factor=0.42,
                orientation="landscape",
            )
            self.assertEqual(fixed_spacing["flipbookReference"]["defaultSpacingMode"], "fixed")
            self.assertEqual(fixed_spacing["flipbookReference"]["defaultSpacingFactor"], 0.42)
            self.assertEqual(fixed_spacing["rotationModel"]["orientation"], "landscape")

    def test_numeric_segid_masks_are_mapped_to_pdb_chain_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            _write_slice(tmp_path / "slice_1_first_frame.pdb", "X", [0.10, 0.20], segid="7")
            (tmp_path / flipbook.MASK_METADATA_FILENAME).write_text(
                "ResidueID,ChainID,Masked\n1,7,True\n2,7,False\n",
                encoding="utf-8",
            )

            manifest = build_molstar_manifest(tmp_path)

            self.assertEqual(manifest["residues"][0]["chain"], "X")
            self.assertEqual(manifest["residues"][0]["segid"], "7")
            self.assertEqual(manifest["maskSummary"]["maskedKeys"], ["X:1"])

    def test_packaged_molstar_assets_are_available(self) -> None:
        self.assertTrue((ASSET_DIR / "script.js").is_file())
        self.assertTrue((ASSET_DIR / "vendor/molstar/5.4.2/molstar.js").is_file())
        self.assertTrue((ASSET_DIR / "vendor/molstar/5.4.2/molstar.css").is_file())
        self.assertTrue((ASSET_DIR / "vendor/molstar/5.4.2/NOTICE.txt").is_file())
        script = (ASSET_DIR / "script.js").read_text(encoding="utf-8")
        self.assertIn('const CAMERA_MODES = new Set(["orthographic", "perspective"])', script)
        self.assertIn('mode: state.cameraMode', script)
        self.assertIn('step="0.01"', script)
        self.assertIn("function roundSpacing(value)", script)
        self.assertIn("spacingStep: spacingStep()", script)
        self.assertIn("function automaticSpacingFactor()", script)
        self.assertIn('defaultSpacingMode === "auto"', script)
        self.assertNotIn("queueGeometryUpdate(true)", script)

    def test_top_level_molstar_helpers_are_lazy_exports(self) -> None:
        from rmsx import build_molstar_manifest as exported_builder

        self.assertIs(exported_builder, build_molstar_manifest)


if __name__ == "__main__":
    unittest.main()
