import tempfile
import unittest
from pathlib import Path
from unittest import mock
from io import StringIO
from contextlib import redirect_stdout


try:
    import numpy as np
    import pandas as pd
    import MDAnalysis as mda

    from rmsx import core as rmsx_core

    HAS_CORE_DEPS = True
except Exception:  # pragma: no cover - lightweight environments skip these tests
    HAS_CORE_DEPS = False


@unittest.skipUnless(HAS_CORE_DEPS, "Masking tests require pandas and MDAnalysis")
class TestMaskingHelpers(unittest.TestCase):
    def test_normalize_mask_selections(self) -> None:
        self.assertEqual(rmsx_core.normalize_mask_selections(None), [])
        self.assertEqual(rmsx_core.normalize_mask_selections(" segid A and resid 1:5 "), ["segid A and resid 1:5"])
        self.assertEqual(
            rmsx_core.normalize_mask_selections(["segid A", "  ", "resid 10:12"]),
            ["segid A", "resid 10:12"],
        )
        with self.assertRaises(TypeError):
            rmsx_core.normalize_mask_selections(123)

    def test_apply_mask_clipping_raises_for_fully_masked_single_chain(self) -> None:
        data_frame = pd.DataFrame(
            {
                "ResidueID": [1, 2],
                "ChainID": ["A", "A"],
                "slice_1.dcd": [1.0, 5.0],
                "slice_2.dcd": [2.0, 6.0],
            }
        )
        mask_metadata = pd.DataFrame(
            {
                "ResidueID": [1, 2],
                "ChainID": ["A", "A"],
                "Masked": [True, True],
            }
        )

        with self.assertRaises(ValueError):
            rmsx_core.apply_mask_clipping(data_frame, mask_metadata)

    def test_summaries_and_global_range_ignore_masked_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            first_dir = tmp_path / "case_one"
            second_dir = tmp_path / "case_two"
            first_dir.mkdir()
            second_dir.mkdir()
            csv_one = first_dir / "one.csv"
            csv_two = second_dir / "two.csv"

            pd.DataFrame(
                {
                    "ResidueID": [1, 2],
                    "ChainID": ["A", "A"],
                    "slice_1.dcd": [1.0, 100.0],
                    "slice_2.dcd": [2.0, 200.0],
                }
            ).to_csv(csv_one, index=False)
            pd.DataFrame(
                {
                    "ResidueID": [3, 4],
                    "ChainID": ["B", "B"],
                    "slice_1.dcd": [3.0, 4.0],
                    "slice_2.dcd": [5.0, 6.0],
                }
            ).to_csv(csv_two, index=False)

            pd.DataFrame(
                {
                    "ResidueID": [1, 2],
                    "ChainID": ["A", "A"],
                    "Masked": [False, True],
                }
            ).to_csv(rmsx_core.get_mask_metadata_path(csv_one), index=False)
            pd.DataFrame(
                {
                    "ResidueID": [3, 4],
                    "ChainID": ["B", "B"],
                    "Masked": [False, False],
                }
            ).to_csv(rmsx_core.get_mask_metadata_path(csv_two), index=False)

            top_df, _ = rmsx_core.summarize_rmsx(csv_one, n=1, print_output=False)
            self.assertEqual(int(top_df.iloc[0]["ResidueID"]), 1)
            self.assertEqual(float(top_df.iloc[0]["RMSX"]), 2.0)

            buffer = StringIO()
            with redirect_stdout(buffer):
                rmsx_core.summarize_rmsx(csv_one, n=1, print_output=True)
            self.assertIn("excluded 1 masked residue(s)", buffer.getvalue())

            global_min, global_max = rmsx_core.compute_global_rmsx_min_max([str(csv_one), str(csv_two)])
            self.assertEqual(global_min, 1.0)
            self.assertEqual(global_max, 6.0)

    def test_create_r_plot_raises_clear_mask_error_on_masked_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rmsx_csv = tmp_path / "rmsx_test_1.000_ns.csv"
            rmsd_csv = tmp_path / "rmsd.csv"
            rmsf_csv = tmp_path / "rmsf.csv"

            pd.DataFrame(
                {
                    "ResidueID": [1, 2],
                    "ChainID": ["A", "A"],
                    "slice_1.dcd": [1.0, 2.0],
                }
            ).to_csv(rmsx_csv, index=False)
            pd.DataFrame({"ResidueID": [1, 2], "ChainID": ["A", "A"], "Masked": [True, False]}).to_csv(
                rmsx_core.get_mask_metadata_path(rmsx_csv),
                index=False,
            )
            pd.DataFrame({"Frame": [0], "Time": [0], "RMSD": [0.0]}).to_csv(rmsd_csv, index=False)
            pd.DataFrame({"ResidueID": [1, 2], "RMSF": [0.0, 0.0]}).to_csv(rmsf_csv, index=False)

            calls = [
                mock.Mock(returncode=0, stdout="", stderr=""),
                mock.Mock(returncode=1, stdout="", stderr="ggpattern unavailable"),
            ]
            with mock.patch.object(rmsx_core.subprocess, "run", side_effect=calls):
                with self.assertRaises(RuntimeError) as ctx:
                    rmsx_core.create_r_plot(
                        str(rmsx_csv),
                        str(rmsd_csv),
                        str(rmsf_csv),
                        rscript_executable="Rscript",
                        verbose=False,
                    )
            self.assertIn("Masked heatmap rendering failed", str(ctx.exception))


@unittest.skipUnless(HAS_CORE_DEPS, "Masking integration tests require pandas and MDAnalysis")
class TestMaskingIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        demo_dir = Path("/Users/finn/Documents/GitHub/rmsx/rmsx/test_files")
        cls.single_topology = demo_dir / "1UBQ.pdb"
        cls.single_trajectory = demo_dir / "mon_sys.dcd"
        cls.multi_topology = demo_dir / "protease_backbone.pdb"
        cls.multi_trajectory = demo_dir / "short_protease_backbone.dcd"

    def test_run_rmsx_and_shift_clip_masked_residues(self) -> None:
        universe = mda.Universe(str(self.single_topology))
        chain_id = str(np.unique(universe.atoms.segids)[0])
        residue_ids = sorted(int(residue.resid) for residue in universe.select_atoms(f"protein and name CA and segid {chain_id}").residues)
        mask_selection = f"segid {chain_id} and resid {residue_ids[0]}:{residue_ids[4]}"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            rmsx_core.run_rmsx(
                topology_file=str(self.single_topology),
                trajectory_file=str(self.single_trajectory),
                output_dir=str(output_dir / "rmsx_case"),
                num_slices=3,
                chain_sele=chain_id,
                verbose=False,
                make_plot=False,
                summary_n=None,
                mask=mask_selection,
            )
            rmsx_csv = next((output_dir / "rmsx_case" / f"chain_{chain_id}_rmsx").glob("rmsx_*.csv"))
            rmsx_df = pd.read_csv(rmsx_csv)
            rmsx_mask = pd.read_csv(rmsx_core.get_mask_metadata_path(rmsx_csv))
            value_columns = [column for column in rmsx_df.columns if column not in ("ResidueID", "ChainID")]
            unmasked_values = rmsx_df.loc[~rmsx_mask["Masked"], value_columns].to_numpy(dtype=float)
            masked_values = rmsx_df.loc[rmsx_mask["Masked"], value_columns].to_numpy(dtype=float)
            self.assertGreater(masked_values.size, 0)
            self.assertGreaterEqual(masked_values.min(), unmasked_values.min())
            self.assertLessEqual(masked_values.max(), unmasked_values.max())

            rmsx_core.run_shift_map(
                topology_file=str(self.single_topology),
                trajectory_file=str(self.single_trajectory),
                output_dir=str(output_dir / "shift_case"),
                num_slices=3,
                chain_sele=chain_id,
                verbose=False,
                make_plot=False,
                summary_n=None,
                mask=mask_selection,
            )
            shift_csv = next((output_dir / "shift_case" / f"chain_{chain_id}_shiftmap").glob("rmsx_*.csv"))
            shift_df = pd.read_csv(shift_csv)
            shift_mask = pd.read_csv(rmsx_core.get_mask_metadata_path(shift_csv))
            shift_columns = [column for column in shift_df.columns if column not in ("ResidueID", "ChainID")]
            unmasked_shift = shift_df.loc[~shift_mask["Masked"], shift_columns].to_numpy(dtype=float)
            masked_shift = shift_df.loc[shift_mask["Masked"], shift_columns].to_numpy(dtype=float)
            self.assertGreater(masked_shift.size, 0)
            self.assertGreaterEqual(masked_shift.min(), unmasked_shift.min())
            self.assertLessEqual(masked_shift.max(), unmasked_shift.max())

    def test_run_rmsx_errors_when_single_chain_is_fully_masked(self) -> None:
        universe = mda.Universe(str(self.single_topology))
        chain_id = str(np.unique(universe.atoms.segids)[0])

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                rmsx_core.run_rmsx(
                    topology_file=str(self.single_topology),
                    trajectory_file=str(self.single_trajectory),
                    output_dir=str(Path(tmpdir) / "fully_masked"),
                    num_slices=2,
                    chain_sele=chain_id,
                    verbose=False,
                    make_plot=False,
                    summary_n=None,
                    mask=f"segid {chain_id}",
                )

    def test_all_chain_rmsx_handles_fully_masked_chain_and_writes_combined_mask_metadata(self) -> None:
        universe = mda.Universe(str(self.multi_topology))
        chain_ids = sorted(str(chain_id) for chain_id in np.unique(universe.atoms.segids))
        self.assertGreaterEqual(len(chain_ids), 2, "Expected packaged protease demo to include multiple chains")
        masked_chain = chain_ids[0]
        unmasked_chain = chain_ids[1]

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(rmsx_core, "create_r_plot", return_value=True):
                combined_dir = rmsx_core.all_chain_rmsx(
                    topology_file=str(self.multi_topology),
                    trajectory_file=str(self.multi_trajectory),
                    output_dir=str(Path(tmpdir) / "multi_case"),
                    num_slices=2,
                    verbose=False,
                    summary_n=None,
                    mask=f"segid {masked_chain}",
                )

            combined_mask = pd.read_csv(Path(combined_dir) / rmsx_core.MASK_METADATA_FILENAME)
            self.assertTrue((combined_mask["Masked"]).any())
            self.assertIn(masked_chain, set(combined_mask["ChainID"].astype(str)))

            masked_csv = next((Path(tmpdir) / "multi_case" / f"chain_{masked_chain}_rmsx").glob("rmsx_*.csv"))
            unmasked_csv = next((Path(tmpdir) / "multi_case" / f"chain_{unmasked_chain}_rmsx").glob("rmsx_*.csv"))

            masked_df = pd.read_csv(masked_csv)
            unmasked_df = pd.read_csv(unmasked_csv)
            masked_meta = pd.read_csv(rmsx_core.get_mask_metadata_path(masked_csv))
            unmasked_meta = pd.read_csv(rmsx_core.get_mask_metadata_path(unmasked_csv))
            value_columns = [column for column in masked_df.columns if column not in ("ResidueID", "ChainID")]

            global_min, global_max = rmsx_core.compute_global_rmsx_min_max([str(masked_csv), str(unmasked_csv)])
            masked_values = masked_df.loc[masked_meta["Masked"], value_columns].to_numpy(dtype=float)
            unmasked_values = unmasked_df.loc[~unmasked_meta["Masked"], value_columns].to_numpy(dtype=float)

            self.assertGreater(masked_values.size, 0)
            self.assertGreaterEqual(masked_values.min(), global_min)
            self.assertLessEqual(masked_values.max(), global_max)
            self.assertGreater(unmasked_values.size, 0)


if __name__ == "__main__":
    unittest.main()
