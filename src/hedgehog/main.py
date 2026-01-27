import re
import subprocess
import warnings
from enum import Enum
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import typer
from rdkit import Chem
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from hedgehog.configs.logger import LoggerSingleton, load_config, logger
from hedgehog.pipeline import calculate_metrics
from hedgehog.utils.data_prep import prepare_input_data
from hedgehog.utils.mol_index import assign_mol_idx

mpl.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")

DEFAULT_CONFIG_PATH = "./src/hedgehog/configs/config.yml"
SAMPLED_MOLS_FILENAME = "sampled_molecules.csv"
STAGE_OVERRIDE_KEY = "_run_single_stage_override"

# Supported SMI-like file extensions for ligand preparation tool
SMI_EXTENSIONS = {"smi", "ismi", "cmi", "txt"}


def _validate_input_path(input_path):
    """Validate input path and return Path object if valid, None otherwise."""
    if "*" in input_path or "?" in input_path:
        return None
    input_path_obj = Path(input_path)
    return input_path_obj if input_path_obj.exists() else None


def _folder_is_empty(folder: Path) -> bool:
    """Check if a folder doesn't exist or has no contents."""
    return not folder.exists() or not any(folder.iterdir())


def _get_unique_results_folder(base_folder) -> Path:
    """
    Generate a unique folder name by appending a number if the folder already exists.

    If the base folder exists and contains any files or directories, creates a new folder
    with an incremented suffix (e.g., results -> results1 -> results2).

    Parameters
    ----------
    base_folder : Path or str
        Path object or string representing the base folder

    Returns
    -------
    Path
        Path object with a unique folder name
    """
    base_folder = Path(base_folder)

    if _folder_is_empty(base_folder):
        return base_folder

    parent = base_folder.parent
    base_name = base_folder.name

    match = re.match(r"^(.+?)(\d+)$", base_name)
    name_prefix = match.group(1) if match else base_name
    counter = int(match.group(2)) + 1 if match else 1

    while True:
        new_folder = parent / f"{name_prefix}{counter}"
        if _folder_is_empty(new_folder):
            return new_folder
        counter += 1


def _canonicalize_smiles(smi: str) -> str | None:
    """Canonicalize a SMILES string using RDKit. Returns None if invalid."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def preprocess_input_with_rdkit(input_path, folder_to_save, log) -> str | None:
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
    log
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
            canonical_smi = _canonicalize_smiles(str(row["smiles"]))
            if canonical_smi is None:
                continue
            row_data = {"smiles": canonical_smi}
            if "model_name" in row:
                row_data["model_name"] = row["model_name"]
            cleaned_data.append(row_data)

        if not cleaned_data:
            return None

        output_df = pd.DataFrame(cleaned_data)
        initial_count = len(output_df)

        output_df = output_df.drop_duplicates(
            subset=["smiles", "model_name"], keep="first"
        ).reset_index(drop=True)

        duplicates_removed = initial_count - len(output_df)
        if duplicates_removed > 0:
            log.info(
                "Removed %d duplicate molecules within models",
                duplicates_removed,
            )

        output_df.to_csv(prepared_output, index=False)
        log.info(
            "RDKit preprocessing: %d molecules saved to %s",
            len(output_df),
            prepared_output,
        )
        return str(prepared_output)
    except Exception as e:
        log.debug("RDKit preprocessing failed: %s", e)
        return None


def _get_input_format_flag(extension: str) -> str | None:
    """Get the input format flag for ligand preparation tool based on file extension."""
    ext = extension.lower().lstrip(".")
    if ext == "csv":
        return "-icsv"
    return "-ismi" if ext in SMI_EXTENSIONS else None


def preprocess_input_with_tool(
    input_path, ligand_preparation_tool, folder_to_save, log
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
    log
        Logger instance

    Returns
    -------
    str or None
        Path to prepared file if preprocessing was done, None otherwise
    """
    input_path_obj = _validate_input_path(input_path)
    if not input_path_obj:
        return None

    input_format = _get_input_format_flag(input_path_obj.suffix)
    if not input_format:
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
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return str(prepared_output) if prepared_output.exists() else None
    except Exception:
        return None


def _display_banner() -> None:
    """Display the HEDGEHOG banner."""
    banner_content = (
        "[bold #B29EEE]🦔 HEDGEHOG[/bold #B29EEE]\n"
        "[dim]Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration[/dim]\n"
        "[dim italic]Developed by Ligand Pro[/dim italic]"
    )
    banner = Panel(banner_content, border_style="#B29EEE", padding=(0, 1), expand=False)
    console.print("")
    console.print(banner)
    console.print("")


def _apply_cli_overrides(
    config_dict: dict,
    generated_mols_path: str | None,
    stage: "Stage | None",
) -> None:
    """Apply CLI argument overrides to config dictionary."""
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


def _resolve_output_folder(
    config_dict: dict,
    reuse_folder: bool,
    force_new_folder: bool,
    stage: "Stage | None",
    generated_mols_path: str | None,
) -> Path:
    """Determine and log the appropriate output folder based on CLI flags."""
    original_folder = Path(config_dict["folder_to_save"])

    if reuse_folder:
        logger.info(
            "[#B29EEE]Folder mode:[/#B29EEE] Reusing folder '%s' (--reuse flag)",
            original_folder,
        )
        return original_folder

    if force_new_folder:
        folder = _get_unique_results_folder(original_folder)
        if folder != original_folder:
            logger.info(
                "[#B29EEE]Folder mode:[/#B29EEE] Creating new folder '%s' "
                "(--force-new flag)",
                folder,
            )
        return folder

    # Auto-mode: reuse for stage reruns, create new otherwise
    if stage and not generated_mols_path:
        logger.info(
            "[#B29EEE]Folder mode:[/#B29EEE] Reusing folder '%s' for stage execution",
            original_folder,
        )
        return original_folder

    folder = _get_unique_results_folder(original_folder)
    if folder != original_folder:
        logger.info(
            "[#B29EEE]Folder mode:[/#B29EEE] Folder '%s' "
            "contains results. Using '%s' instead.",
            original_folder,
            folder,
        )
    return folder


def _preprocess_input(
    config_dict: dict,
    folder_to_save: Path,
) -> None:
    """Preprocess input molecules using ligand preparation tool or RDKit."""
    ligand_preparation_tool = config_dict.get("ligand_preparation_tool")
    original_input_path = config_dict.get("generated_mols_path")

    if not original_input_path:
        return

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


def _save_sampled_molecules(
    data: pd.DataFrame,
    folder_to_save: Path,
    should_save: bool,
) -> None:
    """Save sampled molecules to input directory if requested."""
    if not should_save:
        return

    folder_to_save.mkdir(parents=True, exist_ok=True)
    input_dir = folder_to_save / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_path = input_dir / SAMPLED_MOLS_FILENAME

    data.to_csv(output_path, index=False)
    logger.info(
        "[#B29EEE]✓[/#B29EEE] Sampled total of %d molecules saved to %s",
        len(data),
        output_path,
    )


# Initialize Typer app and Rich console
app = typer.Typer(
    name="HEDGEHOG",
    help=(
        "🦔 Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration - "
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

    @property
    def description(self) -> str:
        """Return human-readable description for this stage."""
        descriptions = {
            Stage.descriptors: "Compute 22 physicochemical descriptors per molecule",
            Stage.struct_filters: "Apply structural filters (Lilly, NIBR, PAINS, etc.)",
            Stage.synthesis: (
                "Evaluate synthetic accessibility using retrosynthesis "
                "(AiZynthFinder) and other metrics"
            ),
            Stage.docking: "Calculate docking scores with Smina/Gnina",
        }
        return descriptions[self]


@app.command()
def run(
    generated_mols_path: str | None = typer.Option(
        None,
        "--mols",
        "-m",
        help="Path or glob pattern to generated SMILES files (overrides config)",
    ),
    stage: Stage | None = typer.Option(
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
    reuse_folder: bool = typer.Option(
        False,
        "--reuse",
        help="Force reuse of existing results folder (useful for reruns)",
    ),
    force_new_folder: bool = typer.Option(
        False,
        "--force-new",
        help="Force creation of a new results folder even when rerunning stages",
    ),
) -> None:
    """
    Run the molecular analysis pipeline.

    Examples
    --------
    \b
    # Run full pipeline (auto-creates new folder if results exist)
    uv run hedgehog run

    \b
    # Alternatively, using the short-name alias:
    uv run hedge run

    \b
    # Rerun specific stage (auto-reuses existing folder)
    uv run hedgehog run --stage docking

    \b
    # Run stage with new molecules (auto-creates new folder)
    uv run hedgehog run --stage descriptors --mols data/*.csv

    \b
    # Force reuse existing folder
    uv run hedgehog run --reuse

    \b
    # Force create new folder even for stage rerun
    uv run hedgehog run --stage docking --force-new
    """
    _display_banner()

    if reuse_folder and force_new_folder:
        logger.error(
            "[red]Error:[/red] Cannot use --reuse and --force-new together. "
            "Please choose one."
        )
        raise typer.Exit(code=1)

    config_dict = load_config(DEFAULT_CONFIG_PATH)
    _apply_cli_overrides(config_dict, generated_mols_path, stage)

    folder_to_save = _resolve_output_folder(
        config_dict, reuse_folder, force_new_folder, stage, generated_mols_path
    )
    config_dict["folder_to_save"] = str(folder_to_save)
    LoggerSingleton().configure_log_directory(folder_to_save)

    _preprocess_input(config_dict, folder_to_save)

    data = prepare_input_data(config_dict, logger)

    if "mol_idx" not in data.columns or data["mol_idx"].isna().all():
        data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)

    should_save = config_dict.get("save_sampled_mols", False) or stage is not None
    _save_sampled_molecules(data, folder_to_save, should_save)

    logger.info("[bold #B29EEE]Starting pipeline...[/bold #B29EEE]")
    calculate_metrics(data, config_dict)
    logger.info(
        "[bold][#B29EEE]Ligand Pro[/#B29EEE] thanks you for using "
        "[#B29EEE]🦔 HEDGEHOG[/#B29EEE]![/bold]"
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

    for stage in Stage:
        table.add_row(stage.value, stage.description)

    console.print(table)
    console.print("\n[dim]Example (1): uv run hedgehog run --stage descriptors[/dim]")
    console.print("[dim]Example (2): uv run hedge run --help [/dim]")


@app.command()
def version() -> None:
    """Display version information."""
    console.print("[bold #B29EEE]🦔 HEDGEHOG[/bold #B29EEE] version [bold]1.0.0[/bold]")
    console.print(
        "[dim]Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration[/dim]"
    )
    console.print("[dim]Developed by [bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim]")


if __name__ == "__main__":
    app()


@app.command()
def tui() -> None:
    """
    Launch the interactive TUI (Text User Interface).

    Requires Node.js to be installed. The TUI provides a visual interface
    for configuring and running the pipeline.

    Examples
    --------
    \b
    uv run hedgehog tui
    """
    import shutil
    import subprocess
    from pathlib import Path

    # Find the TUI directory
    tui_dir = Path(__file__).parent.parent.parent.parent / "tui"

    if not tui_dir.exists():
        console.print(
            "[red]Error:[/red] TUI directory not found. "
            "Please ensure the TUI is installed."
        )
        raise typer.Exit(code=1)

    # Check for Node.js
    node_path = shutil.which("node")
    if not node_path:
        console.print(
            "[red]Error:[/red] Node.js not found. "
            "Please install Node.js to use the TUI."
        )
        raise typer.Exit(code=1)

    # Check if TUI is built
    dist_dir = tui_dir / "dist"
    if not dist_dir.exists():
        console.print("[yellow]TUI not built. Building...[/yellow]")
        try:
            subprocess.run(
                ["npm", "install"],
                cwd=tui_dir,
                check=True,
                capture_output=True,
            )
            subprocess.run(
                ["npm", "run", "build"],
                cwd=tui_dir,
                check=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            console.print(f"[red]Error building TUI:[/red] {e}")
            raise typer.Exit(code=1)

    # Launch the TUI
    try:
        subprocess.run(
            ["node", str(dist_dir / "index.js")],
            cwd=tui_dir,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        console.print(f"[red]TUI exited with error:[/red] {e.returncode}")
        raise typer.Exit(code=e.returncode)
    except KeyboardInterrupt:
        pass
