import tempfile
import unittest
from pathlib import Path
from unittest import mock


try:
    import numpy as np
    import pandas as pd
    import MDAnalysis as mda

    from rmsx import core as rmsx_core
    from rmsx.addons import lddt as lddt_addon

    HAS_LDDT_DEPS = True
except Exception:  # pragma: no cover - lightweight environments skip these tests
    HAS_LDDT_DEPS = False


@unittest.skipUnless(HAS_LDDT_DEPS, "lDDT masking tests require pandas and MDAnalysis")
class TestLddtMasking(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        demo_dir = Path("/Users/finn/Documents/GitHub/rmsx/rmsx/test_files")
        cls.single_topology = demo_dir / "1UBQ.pdb"
        cls.single_trajectory = demo_dir / "mon_sys.dcd"
        cls.multi_topology = demo_dir / "protease_backbone.pdb"
        cls.multi_trajectory = demo_dir / "short_protease_backbone.dcd"

    def _single_chain_mask(self) -> tuple[str, str]:
        universe = mda.Universe(str(self.single_topology))
        chain_id = str(np.unique(universe.atoms.segids)[0])
        residue_ids = sorted(
            int(residue.resid)
            for residue in universe.select_atoms(f"protein and name CA and segid {chain_id}").residues
        )
        return chain_id, f"segid {chain_id} and resid {residue_ids[0]}:{residue_ids[4]}"

    def test_run_lddt_map_clips_masked_values_and_writes_metadata(self) -> None:
        chain_id, mask_selection = self._single_chain_mask()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "lddt_case"

            lddt_addon.run_lddt_map(
                topology_file=str(self.single_topology),
                trajectory_file=str(self.single_trajectory),
                output_dir=str(output_dir),
                num_slices=3,
                chain_sele=chain_id,
                verbose=False,
                interpolate=False,
                make_plot=False,
                summary_n=None,
                mask=mask_selection,
            )

            lddt_csv = next((output_dir / f"chain_{chain_id}_lddtmap").glob("lddt_*.csv"))
            lddt_df = pd.read_csv(lddt_csv)
            mask_df = pd.read_csv(rmsx_core.get_mask_metadata_path(lddt_csv))

            value_columns = [column for column in lddt_df.columns if column not in ("ResidueID", "ChainID")]
            unmasked_values = lddt_df.loc[~mask_df["Masked"], value_columns].to_numpy(dtype=float)
            masked_values = lddt_df.loc[mask_df["Masked"], value_columns].to_numpy(dtype=float)

            self.assertGreater(masked_values.size, 0)
            self.assertTrue((mask_df["Masked"]).any())
            self.assertGreaterEqual(masked_values.min(), unmasked_values.min())
            self.assertLessEqual(masked_values.max(), unmasked_values.max())

    def test_run_lddt_map_errors_when_single_chain_is_fully_masked(self) -> None:
        universe = mda.Universe(str(self.single_topology))
        chain_id = str(np.unique(universe.atoms.segids)[0])

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                lddt_addon.run_lddt_map(
                    topology_file=str(self.single_topology),
                    trajectory_file=str(self.single_trajectory),
                    output_dir=str(Path(tmpdir) / "fully_masked"),
                    num_slices=2,
                    chain_sele=chain_id,
                    verbose=False,
                    interpolate=False,
                    make_plot=False,
                    summary_n=None,
                    mask=f"segid {chain_id}",
                )

    def test_all_chain_lddt_map_uses_unmasked_bounds_and_writes_combined_metadata(self) -> None:
        universe = mda.Universe(str(self.multi_topology))
        chain_ids = sorted(str(chain_id) for chain_id in np.unique(universe.atoms.segids))
        self.assertGreaterEqual(len(chain_ids), 2, "Expected packaged protease demo to include multiple chains")

        masked_chain = chain_ids[0]
        partial_chain = chain_ids[1]
        partial_resids = sorted(
            int(residue.resid)
            for residue in universe.select_atoms(f"protein and name CA and segid {partial_chain}").residues
        )
        partial_mask = f"segid {partial_chain} and resid {partial_resids[0]}:{partial_resids[4]}"

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(lddt_addon, "create_r_plot", return_value=True):
                combined_dir = lddt_addon.all_chain_lddt_map(
                    topology_file=str(self.multi_topology),
                    trajectory_file=str(self.multi_trajectory),
                    output_dir=str(Path(tmpdir) / "multi_case"),
                    num_slices=2,
                    verbose=False,
                    interpolate=False,
                    summary_n=None,
                    sync_color_scale=True,
                    mask=[f"segid {masked_chain}", partial_mask],
                )

            combined_mask = pd.read_csv(Path(combined_dir) / rmsx_core.MASK_METADATA_FILENAME)
            self.assertTrue((combined_mask["Masked"]).any())
            self.assertIn(masked_chain, set(combined_mask["ChainID"].astype(str)))

            masked_csv = next((Path(tmpdir) / "multi_case" / f"chain_{masked_chain}_lddtmap").glob("lddt_*.csv"))
            partial_csv = next((Path(tmpdir) / "multi_case" / f"chain_{partial_chain}_lddtmap").glob("lddt_*.csv"))

            masked_df = pd.read_csv(masked_csv)
            partial_df = pd.read_csv(partial_csv)
            masked_meta = pd.read_csv(rmsx_core.get_mask_metadata_path(masked_csv))
            partial_meta = pd.read_csv(rmsx_core.get_mask_metadata_path(partial_csv))
            value_columns = [column for column in masked_df.columns if column not in ("ResidueID", "ChainID")]

            global_min, global_max = rmsx_core.compute_global_rmsx_min_max([str(masked_csv), str(partial_csv)])
            masked_values = masked_df.loc[masked_meta["Masked"], value_columns].to_numpy(dtype=float)
            partial_unmasked_values = partial_df.loc[~partial_meta["Masked"], value_columns].to_numpy(dtype=float)

            self.assertGreater(masked_values.size, 0)
            self.assertGreater(partial_unmasked_values.size, 0)
            self.assertGreaterEqual(masked_values.min(), global_min)
            self.assertLessEqual(masked_values.max(), global_max)
            self.assertGreaterEqual(masked_values.min(), partial_unmasked_values.min())
            self.assertLessEqual(masked_values.max(), partial_unmasked_values.max())

    def test_run_lddt_flipbook_writes_mask_sidecar_for_run_flipbook(self) -> None:
        universe = mda.Universe(str(self.multi_topology))
        chain_ids = sorted(str(chain_id) for chain_id in np.unique(universe.atoms.segids))
        self.assertGreaterEqual(len(chain_ids), 2, "Expected packaged protease demo to include multiple chains")

        flipbook_calls = []

        def _fake_run_flipbook(directory, **kwargs):
            flipbook_calls.append((directory, kwargs))
            self.assertTrue(
                (Path(directory) / rmsx_core.MASK_METADATA_FILENAME).is_file(),
                "Expected combined masked_residues.csv to be present for FlipBook",
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(lddt_addon, "create_r_plot", return_value=True), \
                 mock.patch.object(lddt_addon, "run_flipbook", side_effect=_fake_run_flipbook):
                combined_dir = lddt_addon.run_lddt_flipbook(
                    topology_file=str(self.multi_topology),
                    trajectory_file=str(self.multi_trajectory),
                    output_dir=str(Path(tmpdir) / "flipbook_case"),
                    num_slices=2,
                    verbose=False,
                    interpolate=False,
                    summary_n=None,
                    sync_color_scale=True,
                    viewer="chimerax",
                    mask=f"segid {chain_ids[0]}",
                )

            self.assertEqual(len(flipbook_calls), 1)
            self.assertEqual(Path(flipbook_calls[0][0]), Path(combined_dir))


if __name__ == "__main__":
    unittest.main()
