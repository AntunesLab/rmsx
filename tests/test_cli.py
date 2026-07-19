from rmsx.cli import main


def test_cli_without_arguments_prints_a_quick_start(capsys) -> None:
    main([])

    output = capsys.readouterr().out
    assert "Quick start:" in output
    assert "rmsx topology.pdb trajectory.dcd" in output
    assert "--no-plot" in output
