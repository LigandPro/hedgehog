import pytest

from hedgehog.docking.configuration import (
    DockingConfigError,
    normalize_docking_config,
    validate_docking_config,
)
from hedgehog.docking.metadata import _parse_tools_config


def test_normalize_resolves_shared_paths_relative_to_config(tmp_path, monkeypatch):
    config_dir = tmp_path / "configs"
    input_dir = config_dir / "inputs"
    input_dir.mkdir(parents=True)
    receptor = input_dir / "receptor.pdb"
    autobox = input_dir / "reference.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    autobox.write_text("$$$$\n", encoding="utf-8")
    monkeypatch.setenv("DOCKING_INPUTS", "inputs")

    normalized = normalize_docking_config(
        {
            "receptor_pdb": "${DOCKING_INPUTS}/receptor.pdb",
            "autobox_ligand": "inputs/reference.sdf",
        },
        config_dir / "config_docking.yml",
        ["gnina"],
    )

    assert normalized["receptor_pdb"] == str(receptor.resolve())
    assert normalized["autobox_ligand"] == str(autobox.resolve())


def test_normalize_ignores_unselected_matcha_environment(tmp_path):
    normalized = normalize_docking_config(
        {
            "receptor_pdb": "receptor.pdb",
            "autobox_ligand": "reference.sdf",
            "matcha_config": {
                "checkpoint_root": "${UNSET_MATCHA_ROOT}",
            },
        },
        tmp_path / "config_docking.yml",
        ["gnina"],
    )

    assert normalized["matcha_config"]["checkpoint_root"] == "${UNSET_MATCHA_ROOT}"


def test_normalize_resolves_matcha_checkpoint_root_without_public_training_config(
    tmp_path,
):
    checkpoint_root = tmp_path / "checkpoints"
    run_dir = checkpoint_root / "model-run"
    run_dir.mkdir(parents=True)
    training_config = run_dir / "config.yaml"
    training_config.write_text("seed: 42\n", encoding="utf-8")
    (tmp_path / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
    (tmp_path / "reference.sdf").write_text("$$$$\n", encoding="utf-8")

    normalized = normalize_docking_config(
        {
            "receptor_pdb": "receptor.pdb",
            "autobox_ligand": "reference.sdf",
            "matcha_config": {
                "backend": "docking",
                "checkpoint_root": str(checkpoint_root),
                "checkpoint_run": "model-run",
            },
        },
        tmp_path / "config_docking.yml",
        ["matcha"],
    )

    assert normalized["matcha_config"]["checkpoint_root"] == str(
        checkpoint_root.resolve()
    )
    assert "training_config" not in normalized["matcha_config"]
    validate_docking_config(
        {
            **normalized,
            "receptor_pdb": str(tmp_path / "receptor.pdb"),
            "autobox_ligand": str(tmp_path / "reference.sdf"),
        },
        ["matcha"],
    )


@pytest.mark.parametrize(
    "option",
    [
        "batch_size",
        "checkpoint_name",
        "concurrency",
        "data_workers",
        "gpus",
        "sample_timeout_seconds",
        "training_config",
        "update_checkout",
        "uv_bin",
    ],
)
def test_validate_rejects_removed_matcha_options(tmp_path, option):
    receptor = tmp_path / "receptor.pdb"
    autobox = tmp_path / "reference.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    autobox.write_text("$$$$\n", encoding="utf-8")

    with pytest.raises(DockingConfigError, match="unsupported matcha_config option"):
        validate_docking_config(
            {
                "receptor_pdb": str(receptor),
                "autobox_ligand": str(autobox),
                "matcha_config": {option: "unused"},
            },
            ["matcha"],
        )


def test_validate_requires_search_box_for_selected_engine(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")

    with pytest.raises(DockingConfigError, match="gnina requires autobox_ligand"):
        validate_docking_config({"receptor_pdb": str(receptor)}, ["gnina"])


def test_validate_accepts_shared_autobox(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    autobox = tmp_path / "reference.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    autobox.write_text("$$$$\n", encoding="utf-8")

    validate_docking_config(
        {
            "receptor_pdb": str(receptor),
            "autobox_ligand": str(autobox),
            "smina_config": {},
            "gnina_config": {},
        },
        ["smina", "gnina"],
    )


def test_writer_inherits_shared_autobox_and_padding(tmp_path):
    from hedgehog.docking.config_writer import _create_per_molecule_configs

    receptor = tmp_path / "receptor.pdb"
    molecule = tmp_path / "molecule.sdf"
    autobox = tmp_path / "reference.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    molecule.write_text("\n", encoding="utf-8")
    autobox.write_text("\n", encoding="utf-8")

    entries = _create_per_molecule_configs(
        {
            "autobox_ligand": str(autobox),
            "autobox_add": 4,
            "gnina_config": {"cpu": 1, "num_modes": 9},
        },
        tmp_path,
        receptor,
        [molecule],
        "gnina",
    )

    config_text = entries[0][1].read_text(encoding="utf-8")
    assert f"autobox_ligand = {autobox}" in config_text
    assert "autobox_add = 4" in config_text
    assert "num_modes = 9" in config_text


def test_writer_rejects_unknown_engine_options(tmp_path):
    from hedgehog.docking.config_writer import _create_per_molecule_configs

    receptor = tmp_path / "receptor.pdb"
    molecule = tmp_path / "molecule.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    molecule.write_text("\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported gnina_config options"):
        _create_per_molecule_configs(
            {"gnina_config": {"typo_exhaustivness": 8}},
            tmp_path,
            receptor,
            [molecule],
            "gnina",
        )


def test_validate_rejects_unknown_top_level_option(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    autobox = tmp_path / "reference.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    autobox.write_text("\n", encoding="utf-8")

    with pytest.raises(DockingConfigError, match="unsupported top-level option"):
        validate_docking_config(
            {
                "receptor_pdb": str(receptor),
                "autobox_ligand": str(autobox),
                "gnina_config": {},
                "matcha_config": {},
                "smina_config": {},
                "typo_parallel_job": 4,
            },
            ["gnina"],
        )


def test_parse_tools_is_explicit_and_deduplicated():
    assert _parse_tools_config({"tools": ["gnina", "smina", "gnina"]}) == [
        "gnina",
        "smina",
    ]


def test_parse_tools_rejects_typos():
    with pytest.raises(ValueError, match="Unsupported docking tool"):
        _parse_tools_config({"tools": "gnnia"})
