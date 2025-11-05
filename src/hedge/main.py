import logging
import subprocess
import warnings
from enum import Enum
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import typer
from rdkit import Chem
from rdkit.Chem import SanitizeMol
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from hedge.configs.logger import load_config, logger
from hedge.pipeline import calculate_metrics
from hedge.utils.data_prep import prepare_input_data
from hedge.utils.mol_index import assign_mol_idx

mpl.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")

DEFAULT_CONFIG_PATH = "./src/hedge/configs/config.yml"
SAMPLED_MOLS_FILENAME = "sampledMols.csv"
STAGE_OVERRIDE_KEY = "_run_single_stage_override"


def _validate_input_path(input_path: str) -> Path | None:
    """Validate input path and return Path object if valid, None otherwise."""
    if "*" in input_path or "?" in input_path:
        return None
    input_path_obj = Path(input_path)
    return input_path_obj if input_path_obj.exists() else None


def preprocess_input_with_rdkit(
    input_path: str,
    folder_to_save: Path,
    logger: logging.Logger,
) -> str | None:
    """
    Preprocess input CSV file using RDKit.

    Fallback when ligand_preparation_tool is not provided.
    Canonicalizes SMILES and removes duplicates within each model.

    Parameters
    ----------
    input_path : str
        Path to the input CSV file
    folder_to_save : Path
        Output folder for prepared file
    logger
        Logger instance

    Returns
    -------
    str or None
        Path to prepared CSV file if preprocessing was done, None otherwise
    """
    input_path_obj = _validate_input_path(input_path)
    if not input_path_obj:
        return None

    try:
        folder_to_save.mkdir(parents=True, exist_ok=True)
        prepared_output = folder_to_save / f"prepared_{input_path_obj.stem}.csv"

        df = pd.read_csv(input_path_obj)
        cleaned_data = []

        for _, row in df.iterrows():
            smi = str(row["smiles"])
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                SanitizeMol(mol)
                canonical_smi = Chem.MolToSmiles(mol)
                row_data = {"smiles": canonical_smi}
                if "model_name" in row:
                    row_data["model_name"] = row["model_name"]
                cleaned_data.append(row_data)
            except Exception:  # noqa: BLE001, S112
                continue

        if not cleaned_data:
            return None

        output_df = pd.DataFrame(cleaned_data)

        initial_count = len(output_df)
        output_df = output_df.drop_duplicates(
            subset=["smiles", "model_name"], keep="first"
        ).reset_index(drop=True)
        duplicates_removed = initial_count - len(output_df)
        if duplicates_removed > 0:
            logger.info(
                "Removed %d duplicate molecules within models",
                duplicates_removed,
            )

        output_df.to_csv(prepared_output, index=False)
        logger.info(
            "RDKit preprocessing: %d molecules saved to %s",
            len(output_df),
            prepared_output,
        )
        return str(prepared_output)
    except Exception as e:  # noqa: BLE001
        logger.debug("RDKit preprocessing failed: %s", e)
        return None


def preprocess_input_with_tool(
    input_path: str,
    ligand_preparation_tool: str,
    folder_to_save: Path,
    logger: logging.Logger,  # noqa: ARG001
) -> str | None:
    """
    Preprocess input file using external ligand preparation tool.

    Parameters
    ----------
    input_path : str
        Path to the input file (CSV or SMI)
    ligand_preparation_tool : str
        Path to the ligand preparation tool
    folder_to_save : Path
        Output folder for prepared file
    logger
        Logger instance

    Returns
    -------
    str or None
        Path to prepared file if preprocessing was done, None otherwise
    """
    input_path_obj = _validate_input_path(input_path)
    if not input_path_obj:
        return None

    input_ext = input_path_obj.suffix.lower().lstrip(".")
    if input_ext == "csv":
        input_format = "-icsv"
    elif input_ext in ("smi", "ismi", "cmi", "txt"):
        input_format = "-ismi"
    else:
        return None

    folder_to_save.mkdir(parents=True, exist_ok=True)
    prepared_output = folder_to_save / f"prepared_{input_path_obj.stem}.csv"
    cmd = [
        ligand_preparation_tool,
        input_format,
        str(input_path),
        "-ucsv",
        str(prepared_output),
    ]

    try:
        subprocess.run(  # noqa: S603
            cmd, capture_output=True, text=True, check=True
        )
        return str(prepared_output) if prepared_output.exists() else None
    except Exception:  # noqa: BLE001
        return None


# Initialize Typer app and Rich console
app = typer.Typer(
    name="HEDGE",
    help=(
        "🦔 Hierarchical Evaluation of Drug GEnerators - "
        "Benchmark pipeline for generative models"
    ),
    add_completion=False,
    rich_markup_mode="rich",
)
console = Console()


class Stage(str, Enum):
    """Available pipeline stages."""

    descriptors = "descriptors"
    struct_filters = "struct_filters"
    synthesis = "synthesis"
    docking = "docking"


@app.command()
def run(
    generated_mols_path: str | None = typer.Option(
        None,
        "--mols",
        "-m",
        help="Path or glob pattern to generated SMILES files (overrides config)",
    ),
    stage: Stage | None = typer.Option(  # noqa: B008
        None,
        "--stage",
        "-s",
        help=(
            "Run only a specific pipeline stage. Please provide also --mols "
            "argument to specify the molecules path. If no --mols provided, "
            "the pipeline will use the molecules from the config file."
        ),
        case_sensitive=False,
    ),
) -> None:
    r"""
    Run the molecular analysis pipeline.

    Examples
    --------
    \b
    # Run full pipeline
    uv run python main.py run

    \b
    # Run only descriptors stage
    uv run python main.py run --stage descriptors

    \b
    # Override molecules path
    uv run python main.py run --mols data/*.csv
    """
    # Display banner
    console.print(
        Panel.fit(
            "[bold #B29EEE]🦔 HEDGE[/bold #B29EEE]\n"
            "[dim]Hierarchical Evaluation of Drug GEnerators[/dim]\n"
            "[dim italic]Developed by "
            "[bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim italic]",
            border_style="#B29EEE",
        )
    )

    config_dict = load_config(DEFAULT_CONFIG_PATH)

    if generated_mols_path:
        config_dict["generated_mols_path"] = generated_mols_path
        logger.info(
            "[#B29EEE]Override:[/#B29EEE] Using molecules from: %s",
            generated_mols_path,
        )
    elif stage:
        logger.info(
            "Using molecules from config: %s",
            config_dict["generated_mols_path"],
        )

    if stage:
        config_dict[STAGE_OVERRIDE_KEY] = stage.value
        logger.info(
            "[#B29EEE]Override:[/#B29EEE] Running only stage: [bold]%s[/bold]",
            stage.value,
        )

    folder_to_save = Path(config_dict["folder_to_save"])

    ligand_preparation_tool = config_dict.get("ligand_preparation_tool")
    original_input_path = config_dict.get("generated_mols_path") or generated_mols_path
    if original_input_path:
        if ligand_preparation_tool:
            prepared_path = preprocess_input_with_tool(
                original_input_path,
                ligand_preparation_tool,
                folder_to_save,
                logger,
            )
        else:
            prepared_path = preprocess_input_with_rdkit(
                original_input_path, folder_to_save, logger
            )

        if prepared_path:
            config_dict["generated_mols_path"] = prepared_path
            logger.info("Using preprocessed input: %s", prepared_path)

    data = prepare_input_data(config_dict, logger)

    if "mol_idx" not in data.columns or data["mol_idx"].isna().all():
        try:
            data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)
        except Exception as e:
            logger.error("Failed to assign mol_idx: %s", e)
            raise

    save_mols = config_dict.get("save_sampled_mols", False) or stage is not None
    if save_mols:
        folder_to_save.mkdir(parents=True, exist_ok=True)
        output_path = folder_to_save / SAMPLED_MOLS_FILENAME
        data.to_csv(output_path, index=False)
        logger.info(
            "[#B29EEE]✓[/#B29EEE] Sampled total of %d molecules saved to %s",
            len(data),
            output_path,
        )

    # Run metrics calculation pipeline
    logger.info("[bold #B29EEE]Starting pipeline...[/bold #B29EEE]")
    calculate_metrics(data, config_dict)
    logger.info(
        "[bold][#B29EEE]Ligand Pro[/#B29EEE] thanks you for using "
        "[#B29EEE]🦔 HEDGE[/#B29EEE]![/bold]"
    )


@app.command()
def info() -> None:
    """Display pipeline information and available stages."""
    table = Table(
        title="Available Pipeline Stages",
        show_header=True,
        header_style="bold #B29EEE",
    )
    table.add_column("Stage", style="#B29EEE", no_wrap=True)
    table.add_column("Description", style="white")

    stages_info = {
        "descriptors": "Compute 22 physicochemical descriptors per molecule",
        "struct_filters": ("Apply structural filters (Lilly, NIBR, PAINS, etc.)"),
        "synthesis": (
            "Evaluate synthetic accessibility using retrosynthesis "
            "(AiZynthFinder) and other metrics"
        ),
        "docking": "Calculate docking scores with Smina/Gnina",
    }

    for stage_name, description in stages_info.items():
        table.add_row(stage_name, description)

    console.print(table)
    console.print("\n[dim]Example (1): uv run hedge run --stage descriptors[/dim]")
    console.print("[dim]Example (2): uv run hedge run --help [/dim]")


@app.command()
def version() -> None:
    """Display version information."""
    console.print("[bold #B29EEE]🦔 HEDGE[/bold #B29EEE] version [bold]1.0.0[/bold]")
    console.print("[dim]Hierarchical Evaluation of Drug GEnerators[/dim]")
    console.print("[dim]Developed by [bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim]")


if __name__ == "__main__":
    app()
