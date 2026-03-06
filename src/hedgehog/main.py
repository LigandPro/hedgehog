import os
import subprocess
import time
import warnings
from enum import Enum
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import typer
from rdkit import Chem
from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from hedgehog._constants import (
    CFG_DESCRIPTORS,
    CFG_DOCKING,
    CFG_STRUCT_FILTERS,
    KEY_FOLDER_TO_SAVE,
)
from hedgehog.configs.logger import LoggerSingleton, load_config, logger
from hedgehog.pipeline import calculate_metrics
from hedgehog.utils.data_prep import prepare_input_data
from hedgehog.utils.mol_index import assign_mol_idx

mpl.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")

DEFAULT_CONFIG_PATH = str(Path(__file__).resolve().parent / "configs" / "config.yml")
SAMPLED_MOLS_FILENAME = "sampled_molecules.csv"
STAGE_OVERRIDE_KEY = "_run_single_stage_override"
PLAIN_OUTPUT_ENV = "HEDGEHOG_PLAIN_OUTPUT"

# Supported SMI-like file extensions for ligand preparation tool
SMI_EXTENSIONS = {"smi", "ismi", "cmi", "txt"}


def _plain_output_enabled() -> bool:
    return os.environ.get(PLAIN_OUTPUT_ENV, "").strip() == "1"


def _build_progress_columns(console_width: int):
    """Build adaptive progress columns for narrow terminals/tmux panes."""
    if console_width < 120:
        return [
            SpinnerColumn(style="dim"),
            TextColumn("[bold]{task.description}[/bold]"),
            BarColumn(bar_width=20),
            TextColumn("{task.fields[done_total]}"),
            TimeElapsedColumn(),
        ]

    return [
        SpinnerColumn(style="dim"),
        TextColumn("[bold]{task.description}[/bold]"),
        BarColumn(bar_width=40),
        TextColumn("{task.fields[done_total]}"),
        TextColumn("rate {task.fields[rate]}"),
        TextColumn("eta {task.fields[eta]}"),
        TimeElapsedColumn(),
    ]


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
    Generate a unique folder name with sequential number suffix.

    Creates folders like: results/run_1, results/run_2, etc.

    Parameters
    ----------
    base_folder : Path or str
        Path object or string representing the base folder (e.g., results/run)

    Returns
    -------
    Path
        Path object with a unique sequentially numbered folder name
    """
    import re

    base_folder = Path(base_folder)
    if _folder_is_empty(base_folder):
        return base_folder
    parent = base_folder.parent
    base_name = base_folder.name

    # Find existing folders matching pattern {base_name}_N
    max_number = 0
    pattern = re.compile(rf"^{re.escape(base_name)}_(\d+)$")

    if parent.exists():
        for item in parent.iterdir():
            if item.is_dir():
                match = pattern.match(item.name)
                if match:
                    number = int(match.group(1))
                    max_number = max(max_number, number)

    new_folder = parent / f"{base_name}_{max_number + 1}"
    return new_folder


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
    For CSV inputs, performs lightweight schema normalization and duplicate
    removal. Full molecule standardization is handled by the MolPrep stage.

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
        lower_cols = {c.lower(): c for c in df.columns}
        smiles_col = lower_cols.get("smiles")
        if not smiles_col:
            log.debug("CSV preprocessing: missing SMILES column in %s", input_path)
            return None
        if smiles_col != "smiles":
            df = df.rename(columns={smiles_col: "smiles"})

        model_col = lower_cols.get("model_name") or lower_cols.get("name")
        if model_col:
            if model_col != "model_name":
                df = df.rename(columns={model_col: "model_name"})
        else:
            df["model_name"] = input_path_obj.stem

        output_df = df[["smiles", "model_name"]].dropna(subset=["smiles"]).copy()
        initial_count = len(output_df)
        output_df = output_df.drop_duplicates(
            subset=["smiles", "model_name"]
        ).reset_index(drop=True)

        duplicates_removed = initial_count - len(output_df)
        if duplicates_removed > 0:
            log.info(
                "Removed %d duplicate molecules within models",
                duplicates_removed,
            )

        output_df.to_csv(prepared_output, index=False)
        log.info(
            "CSV preprocessing: %d molecules saved to %s",
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
        "-ocsv",
        str(prepared_output),
        "-WAIT",
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        if not prepared_output.exists():
            return None

        # Normalize the preparation tool output to the minimal schema expected by
        # HEDGEHOG input loading: (smiles, model_name).
        df = pd.read_csv(prepared_output)
        lower_cols = {c.lower(): c for c in df.columns}
        smiles_col = lower_cols.get("smiles")
        if not smiles_col:
            log.debug(
                "Ligand preparation output missing SMILES column in %s (cols=%s)",
                prepared_output,
                list(df.columns),
            )
            return None

        if smiles_col != "smiles":
            df = df.rename(columns={smiles_col: "smiles"})

        # Some external preparation tools may emit a per-molecule name/title
        # column; do not treat it as a model identifier. Pin model_name to the
        # input file stem.
        df["model_name"] = input_path_obj.stem
        df = df[["smiles", "model_name"]].dropna(subset=["smiles"])
        df.to_csv(prepared_output, index=False)

        return str(prepared_output)
    except Exception as e:
        log.debug("Ligand preparation failed: %s", e)
        return None


def _display_banner() -> None:
    """Display the HEDGEHOG banner."""
    if _plain_output_enabled():
        return
    banner_content = (
        "[bold]🦔 HEDGEHOG[/bold]\n"
        "[dim]Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration[/dim]\n"
        "[dim italic]Developed by Ligand Pro[/dim italic]"
    )
    banner = Panel(banner_content, border_style="dim", padding=(0, 1), expand=False)
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
            "[bold]Override:[/bold] Using molecules from: %s",
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
            "[bold]Override:[/bold] Running only stage: [bold]%s[/bold]",
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
    original_folder = Path(config_dict[KEY_FOLDER_TO_SAVE])

    if reuse_folder:
        logger.info(
            "[bold]Folder mode:[/bold] Reusing folder '%s' (--reuse flag)",
            original_folder,
        )
        return original_folder

    if force_new_folder:
        folder = _get_unique_results_folder(original_folder)
        if folder != original_folder:
            logger.info(
                "[bold]Folder mode:[/bold] Creating new folder '%s' (--force-new flag)",
                folder,
            )
        return folder

    # Auto-mode: reuse for stage reruns, create new otherwise
    if stage and not generated_mols_path:
        logger.info(
            "[bold]Folder mode:[/bold] Reusing folder '%s' for stage execution",
            original_folder,
        )
        return original_folder

    folder = _get_unique_results_folder(original_folder)
    if folder != original_folder:
        logger.info(
            "[bold]Folder mode:[/bold] Folder '%s' "
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

    input_path_obj = _validate_input_path(original_input_path)
    # If this is a glob pattern or doesn't exist, defer to prepare_input_data().
    if not input_path_obj:
        return

    # Prefer the lightweight RDKit preprocessing for CSV inputs. The external
    # ligand preparation tool is primarily useful for SMI-like inputs.
    if input_path_obj.suffix.lower() == ".csv":
        prepared_path = preprocess_input_with_rdkit(
            original_input_path, folder_to_save, logger
        )
    elif ligand_preparation_tool:
        prepared_path = preprocess_input_with_tool(
            original_input_path,
            ligand_preparation_tool,
            folder_to_save,
            logger,
        )
    else:
        prepared_path = None

    if prepared_path:
        config_dict["generated_mols_path"] = prepared_path
        logger.info("Using preprocessed input: %s", prepared_path)


class CliProgressTracker:
    """Manages Rich progress display for CLI pipeline execution.

    Wraps a :class:`rich.progress.Progress` bar and translates pipeline
    progress events (dicts with ``type``, ``stage``, ``current``, ``total``,
    etc.) into visual updates.  Intended to be used as a context manager so
    that the underlying ``Progress`` widget is properly started/stopped.
    """

    _SHORT_STAGE_NAMES: dict[str, str] = {
        "mol_prep": "Prep",
        "descriptors": "Descriptors",
        "struct_filters": "StructFilters",
        "synthesis": "Synthesis",
        "docking": "Docking",
        "docking_filters": "DockFilters",
        "final_descriptors": "FinalDesc",
        "report": "Report",
    }

    def __init__(self, console: Console) -> None:
        columns = _build_progress_columns(console.width)
        self._progress = Progress(
            *columns,
            refresh_per_second=4,
            console=console,
        )
        self._active_task_id: int | None = None
        self._stage_started_at: float | None = None
        self._current_stage_done: int | None = None
        self._current_stage_total: int | None = None

    # -- context manager --------------------------------------------------

    def __enter__(self) -> "CliProgressTracker":
        self._progress.__enter__()
        return self

    def __exit__(self, *exc: object) -> None:
        self._progress.__exit__(*exc)

    # -- public API -------------------------------------------------------

    def handle_event(self, event: dict) -> None:
        """Handle a progress event emitted by the pipeline."""
        event_type = event.get("type")
        stage = str(event.get("stage", ""))
        stage_index = int(event.get("stage_index", 0) or 0)
        total_stages = int(event.get("total_stages", 0) or 0)
        short_name = self._short_stage_name(stage)
        message = event.get("message")
        task_id = self._get_or_create_task(
            stage_index, total_stages, short_name, message
        )
        if task_id is None:
            return

        if event_type == "stage_start":
            self._handle_stage_start(
                task_id,
                stage_index,
                total_stages,
                short_name,
                message,
                molecules_in=event.get("molecules_in"),
            )
        elif event_type == "stage_progress":
            self._handle_stage_progress(
                event, task_id, stage_index, total_stages, short_name, message
            )
        elif event_type == "stage_complete":
            self._handle_stage_complete(
                task_id,
                stage_index,
                total_stages,
                short_name,
                message,
                molecules_in=event.get("molecules_in"),
                molecules_out=event.get("molecules_out"),
            )

    # -- private helpers --------------------------------------------------

    def _short_stage_name(self, stage: str) -> str:
        if stage in self._SHORT_STAGE_NAMES:
            return self._SHORT_STAGE_NAMES[stage]
        return stage.replace("_", " ").title()

    @staticmethod
    def _format_seconds(value: float | None) -> str:
        if value is None:
            return "-"
        if value < 60:
            return f"{value:.1f}s"
        minutes, seconds = divmod(int(round(value)), 60)
        if minutes < 60:
            return f"{minutes:02d}:{seconds:02d}"
        hours, minutes = divmod(minutes, 60)
        return f"{hours:d}:{minutes:02d}:{seconds:02d}"

    @staticmethod
    def _format_count(value: int | None) -> str:
        if value is None:
            return "-"
        return f"{value:,}"

    @staticmethod
    def _short_progress_message(value: str | None) -> str:
        if not value:
            return ""
        text = " ".join(str(value).split())
        if len(text) > 36:
            return text[:33] + "..."
        return text

    @staticmethod
    def _strip_stage_prefix(message: str, short_name: str) -> str:
        lowered = message.lower()
        stage_lower = short_name.lower()
        if lowered == stage_lower:
            return ""
        for separator in (":", "-", "·", " "):
            prefix = f"{stage_lower}{separator}"
            if lowered.startswith(prefix):
                return message[len(short_name + separator) :].strip()
        return message

    def _progress_description(
        self,
        stage_index: int,
        total_stages: int,
        short_name: str,
        message: str | None = None,
    ) -> str:
        base = f"{stage_index}/{total_stages} - {short_name}"
        short_message = self._short_progress_message(message)
        if short_message:
            short_message = self._strip_stage_prefix(short_message, short_name)
        if short_message:
            return f"{base} · {short_message}"
        return base

    def _get_or_create_task(
        self,
        stage_index: int,
        total_stages: int,
        short_name: str,
        message: str | None = None,
    ) -> int | None:
        if stage_index <= 0:
            return None
        if self._active_task_id is None:
            self._active_task_id = self._progress.add_task(
                self._progress_description(
                    stage_index, total_stages, short_name, message
                ),
                total=100,
                done_total="-/-",
                rate="-",
                eta="-",
            )
        return self._active_task_id

    def _handle_stage_start(
        self,
        task_id: int,
        stage_index: int,
        total_stages: int,
        short_name: str,
        message: str | None,
        molecules_in: int | None = None,
    ) -> None:
        self._stage_started_at = time.perf_counter()
        input_total = int(molecules_in) if molecules_in is not None else None
        if input_total is not None and input_total >= 0:
            self._current_stage_done = 0
            self._current_stage_total = input_total
            done_total = f"0/{self._format_count(input_total)}"
        else:
            self._current_stage_done = None
            self._current_stage_total = None
            done_total = "-/-"
        self._progress.update(
            task_id,
            completed=0,
            description=self._progress_description(
                stage_index, total_stages, short_name, message
            ),
            done_total=done_total,
            rate="-",
            eta="-",
        )

    def _handle_stage_progress(
        self,
        event: dict,
        task_id: int,
        stage_index: int,
        total_stages: int,
        short_name: str,
        message: str | None,
    ) -> None:
        current = int(event.get("current", 0) or 0)
        total = int(event.get("total", 0) or 0)
        self._current_stage_done = current if current >= 0 else None
        self._current_stage_total = total if total > 0 else None
        pct = int(round((current / total) * 100)) if total > 0 else 0
        pct = max(0, min(100, pct))
        elapsed_seconds = (
            time.perf_counter() - self._stage_started_at
            if self._stage_started_at is not None
            else None
        )
        left_count = max(total - current, 0) if total > 0 and current >= 0 else None
        rate = (
            (current / elapsed_seconds)
            if elapsed_seconds and elapsed_seconds > 0 and current > 0
            else None
        )
        eta_seconds = (
            (left_count / rate)
            if left_count is not None and rate and rate > 0
            else None
        )
        self._progress.update(
            task_id,
            completed=pct,
            description=self._progress_description(
                stage_index, total_stages, short_name, message
            ),
            done_total=(
                f"{self._format_count(current)}/{self._format_count(total)}"
                if total > 0
                else "-/-"
            ),
            rate=(f"{rate:,.0f}/s" if rate is not None else "-"),
            eta=self._format_seconds(eta_seconds),
        )

    def _handle_stage_complete(
        self,
        task_id: int,
        stage_index: int,
        total_stages: int,
        short_name: str,
        message: str | None,
        molecules_in: int | None = None,
        molecules_out: int | None = None,
    ) -> None:
        done_total = "-/-"
        stage_total = self._current_stage_total
        stage_done = self._current_stage_done
        if stage_total is not None:
            done = stage_done
            if done is None or done < 0:
                done = stage_total
            done_total = f"{self._format_count(done)}/{self._format_count(stage_total)}"
        else:
            in_count = int(molecules_in) if molecules_in is not None else None
            out_count = int(molecules_out) if molecules_out is not None else None
            if in_count is not None and in_count >= 0:
                if out_count is None or out_count < 0:
                    out_count = in_count
                done_total = (
                    f"{self._format_count(out_count)}/{self._format_count(in_count)}"
                )
            elif out_count is not None and out_count >= 0:
                done_total = (
                    f"{self._format_count(out_count)}/{self._format_count(out_count)}"
                )
        self._progress.update(
            task_id,
            completed=100,
            description=self._progress_description(
                stage_index, total_stages, short_name, message
            ),
            done_total=done_total,
            rate="-",
            eta="-",
        )


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
        "Sampled total of %d molecules saved to %s",
        len(data),
        output_path,
    )


# Initialize Typer app and Rich console
app = typer.Typer(
    name="HEDGEHOG",
    help=(
        "Run the molecular analysis pipeline.\n\n"
        "Examples:\n"
        "  uv run hedgehog\n"
        "  uv run hedge\n"
        "  uv run hedgehog --stage docking\n"
        '  uv run hedgehog --stage descriptors --mols "input/*.csv"\n'
        "  uv run hedgehog --reuse\n"
        "  uv run hedgehog --stage docking --force-new\n"
        "  uv run hedgehog --auto-install\n"
        "  uv run hedgehog --progress\n\n"
        "Stage keys:\n"
        "  mol_prep           Standardize input molecules\n"
        "  descriptors        Compute molecular descriptors\n"
        "  struct_filters     Apply structural alert filters\n"
        "  synthesis          Evaluate synthetic accessibility\n"
        "  docking            Run docking engines\n"
        "  docking_filters    Filter docking poses and interactions\n"
        "  final_descriptors  Recompute descriptors for final set"
    ),
    add_completion=False,
    rich_markup_mode="rich",
)

setup_app = typer.Typer(help="Install optional external tools and assets.")
app.add_typer(setup_app, name="setup")
console = Console(no_color=_plain_output_enabled())


class Stage(str, Enum):
    """Available pipeline stages."""

    mol_prep = "mol_prep"
    descriptors = "descriptors"
    struct_filters = "struct_filters"
    synthesis = "synthesis"
    docking = "docking"
    docking_filters = "docking_filters"
    final_descriptors = "final_descriptors"

    @property
    def description(self) -> str:
        """Return human-readable description for this stage."""
        descriptions = {
            Stage.mol_prep: (
                "Standardize molecules (salts/fragments removal, uncharge, "
                "tautomer canonicalization, stereo removal)"
            ),
            Stage.descriptors: "Compute 22 physicochemical descriptors per molecule",
            Stage.struct_filters: "Apply structural filters (Lilly, NIBR, PAINS, etc.)",
            Stage.synthesis: (
                "Evaluate synthetic accessibility using retrosynthesis "
                "(AiZynthFinder) and other metrics"
            ),
            Stage.docking: "Calculate docking scores with Smina/Gnina",
            Stage.docking_filters: "Filter docking poses by quality and interactions",
            Stage.final_descriptors: "Recompute descriptors on the final filtered set",
        }
        return descriptions[self]


def _run_pipeline_command(
    *,
    config_path: str,
    generated_mols_path: str | None,
    out_dir: str | None,
    stage: Stage | None,
    reuse_folder: bool,
    force_new_folder: bool,
    auto_install: bool,
    show_progress: bool,
) -> None:
    _display_banner()

    if auto_install:
        os.environ["HEDGEHOG_AUTO_INSTALL"] = "1"

    if reuse_folder and force_new_folder:
        logger.error(
            "[red]Error:[/red] Cannot use --reuse and --force-new together. "
            "Please choose one."
        )
        raise typer.Exit(code=1)

    if out_dir and (reuse_folder or force_new_folder):
        logger.error(
            "[red]Error:[/red] Cannot use --out together with --reuse/--force-new. "
            "Please manage uniqueness via --out."
        )
        raise typer.Exit(code=1)

    config_dict = load_config(config_path)
    _apply_cli_overrides(config_dict, generated_mols_path, stage)

    if out_dir:
        folder_to_save = Path(out_dir).resolve()
        logger.info(
            "[bold]Folder override:[/bold] Using output folder '%s' (--out)",
            folder_to_save,
        )
    else:
        folder_to_save = _resolve_output_folder(
            config_dict, reuse_folder, force_new_folder, stage, generated_mols_path
        )
    config_dict[KEY_FOLDER_TO_SAVE] = str(folder_to_save)
    LoggerSingleton().configure_log_directory(folder_to_save)

    _preprocess_input(config_dict, folder_to_save)

    data = prepare_input_data(config_dict, logger)

    if "mol_idx" not in data.columns or data["mol_idx"].isna().all():
        data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)

    should_save = config_dict.get("save_sampled_mols", False) or stage is not None
    _save_sampled_molecules(data, folder_to_save, should_save)

    logger.info("[bold]Starting pipeline...[/bold]")

    shared_console = LoggerSingleton().console

    if _plain_output_enabled() or not show_progress:
        success = calculate_metrics(data, config_dict, None)
    else:
        with CliProgressTracker(shared_console) as tracker:
            success = calculate_metrics(data, config_dict, tracker.handle_event)

    if not success:
        logger.error("Pipeline completed with failures")
        raise typer.Exit(code=1)
    logger.info("Ligand Pro thanks you for using HEDGEHOG!")


@app.callback(invoke_without_command=True)
def run(
    ctx: typer.Context = None,
    config_path: str = typer.Option(
        DEFAULT_CONFIG_PATH,
        "--config",
        "-c",
        help="Master YAML config path (default: src/hedgehog/configs/config.yml).",
        show_default=False,
    ),
    generated_mols_path: str | None = typer.Option(
        None,
        "--mols",
        "-m",
        help="SMILES file path or glob (overrides config).",
    ),
    out_dir: str | None = typer.Option(
        None,
        "--out",
        "-o",
        help="Output directory (overrides config folder_to_save).",
    ),
    stage: Stage | None = typer.Option(
        None,
        "--stage",
        "-s",
        help="Run one stage only; see `hedgehog info`.",
        case_sensitive=False,
        metavar="STAGE",
        show_choices=False,
    ),
    reuse_folder: bool = typer.Option(
        False,
        "--reuse",
        help="Reuse existing results folder.",
    ),
    force_new_folder: bool = typer.Option(
        False,
        "--force-new",
        help="Always create a new results folder.",
    ),
    auto_install: bool = typer.Option(
        False,
        "--auto-install",
        help="Install missing optional tools automatically.",
    ),
    show_progress: bool = typer.Option(
        False,
        "--progress",
        help="Show live progress bar.",
    ),
) -> None:
    """
    Run the molecular analysis pipeline.

    Examples:
    \b
      uv run hedgehog

    \b
      uv run hedge

    \b
      uv run hedgehog --stage docking

    \b
      uv run hedgehog --stage descriptors --mols "input/*.csv"

    \b
      uv run hedgehog --reuse

    \b
      uv run hedgehog --stage docking --force-new

    \b
      uv run hedgehog --auto-install

    \b
      uv run hedgehog --progress

    """
    if ctx is None or ctx.invoked_subcommand is None:
        _run_pipeline_command(
            config_path=config_path,
            generated_mols_path=generated_mols_path,
            out_dir=out_dir,
            stage=stage,
            reuse_folder=reuse_folder,
            force_new_folder=force_new_folder,
            auto_install=auto_install,
            show_progress=show_progress,
        )


@app.command("run")
def run_command(
    config_path: str = typer.Option(
        DEFAULT_CONFIG_PATH,
        "--config",
        "-c",
        help="Master YAML config path (default: src/hedgehog/configs/config.yml).",
        show_default=False,
    ),
    generated_mols_path: str | None = typer.Option(
        None,
        "--mols",
        "-m",
        help="SMILES file path or glob (overrides config).",
    ),
    out_dir: str | None = typer.Option(
        None,
        "--out",
        "-o",
        help="Output directory (overrides config folder_to_save).",
    ),
    stage: Stage | None = typer.Option(
        None,
        "--stage",
        "-s",
        help="Run one stage only; see `hedgehog info`.",
        case_sensitive=False,
        metavar="STAGE",
        show_choices=False,
    ),
    reuse_folder: bool = typer.Option(
        False,
        "--reuse",
        help="Reuse existing results folder.",
    ),
    force_new_folder: bool = typer.Option(
        False,
        "--force-new",
        help="Always create a new results folder.",
    ),
    auto_install: bool = typer.Option(
        False,
        "--auto-install",
        help="Install missing optional tools automatically.",
    ),
    show_progress: bool = typer.Option(
        False,
        "--progress",
        help="Show live progress bar.",
    ),
) -> None:
    """Run the molecular analysis pipeline as an explicit subcommand."""
    _run_pipeline_command(
        config_path=config_path,
        generated_mols_path=generated_mols_path,
        out_dir=out_dir,
        stage=stage,
        reuse_folder=reuse_folder,
        force_new_folder=force_new_folder,
        auto_install=auto_install,
        show_progress=show_progress,
    )


@app.command()
def report(
    results_dir: str = typer.Argument(
        ...,
        help="Path to existing pipeline results directory (e.g., results/run_10)",
    ),
) -> None:
    """
    Regenerate the HTML report from an existing pipeline run directory.

    This command re-generates the report without re-running any pipeline stages.
    Useful when report templates or plotting logic have been updated.

    Examples
    --------
    \b
    uv run hedgehog report results/run_10
    """
    import yaml

    from hedgehog.reporting.report_generator import ReportGenerator

    _display_banner()

    results_path = Path(results_dir).resolve()
    if not results_path.exists():
        console.print(
            f"[red]Error:[/red] Results directory does not exist: {results_path}"
        )
        raise typer.Exit(code=1)

    stages_path = results_path / "stages"
    if not stages_path.exists():
        console.print(
            f"[red]Error:[/red] No stages/ subdirectory found in {results_path}"
        )
        raise typer.Exit(code=1)

    # Load config from resolved master config
    config_path = results_path / "configs" / "master_config_resolved.yml"
    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    config[KEY_FOLDER_TO_SAVE] = str(results_path)

    # Point sub-config keys to actual files in the run's configs directory
    sub_config_keys = [
        CFG_DESCRIPTORS,
        CFG_STRUCT_FILTERS,
        "config_synthesis",
        CFG_DOCKING,
        "config_docking_filters",
        "config_moleval",
    ]
    configs_dir = results_path / "configs"
    for key in sub_config_keys:
        candidate = configs_dir / f"{key}.yml"
        if candidate.exists():
            config[key] = str(candidate)

    # Determine initial molecule count
    initial_count = 0
    sampled_path = results_path / "input" / "sampled_molecules.csv"
    if sampled_path.exists():
        try:
            initial_count = len(pd.read_csv(sampled_path))
        except Exception:
            pass
    if initial_count == 0:
        # Try first stage input CSV as fallback
        for stage_dir in sorted(stages_path.iterdir()):
            if stage_dir.is_dir():
                for csv_file in stage_dir.glob("*.csv"):
                    try:
                        initial_count = len(pd.read_csv(csv_file))
                        break
                    except Exception:
                        continue
            if initial_count > 0:
                break

    # Determine final molecule count
    final_count = 0
    final_path = results_path / "output" / "final_molecules.csv"
    if final_path.exists():
        try:
            final_count = len(pd.read_csv(final_path))
        except Exception:
            pass

    logger.info(
        "Regenerating report for [bold]%s[/bold] (initial=%d, final=%d molecules)",
        results_path,
        initial_count,
        final_count,
    )

    try:
        report_gen = ReportGenerator(
            base_path=results_path,
            stages=[],
            config=config,
            initial_count=initial_count,
            final_count=final_count,
        )
        report_path = report_gen.generate()
        logger.info("[bold]Report generated:[/bold] %s", report_path)
        console.print(f"\n[bold]Report saved to:[/bold] {report_path}\n")
    except Exception as e:
        console.print(f"[red]Error generating report:[/red] {e}")
        raise typer.Exit(code=1) from e


@app.command()
def info() -> None:
    """Display pipeline information and available stages."""
    table = Table(
        title="Available Pipeline Stages",
        show_header=True,
        header_style="bold",
    )
    table.add_column("Stage", style="bold", no_wrap=True)
    table.add_column("Description", style="white")

    for stage in Stage:
        table.add_row(stage.value, stage.description)

    console.print(table)
    console.print("\n[dim]Example (1): uv run hedgehog --stage descriptors[/dim]")
    console.print("[dim]Example (2): uv run hedge --help [/dim]")


@app.command()
def version() -> None:
    """Display version information."""
    if _plain_output_enabled():
        console.print("HEDGEHOG version 1.1.11")
    else:
        console.print("[bold]🦔 HEDGEHOG[/bold] version [bold]1.1.11[/bold]")
    console.print(
        "[dim]Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration[/dim]"
    )
    console.print("[dim]Developed by [bold]Ligand Pro[/bold][/dim]")


@setup_app.command("aizynthfinder")
def setup_aizynthfinder(
    yes: bool = typer.Option(
        True,
        "--yes/--no-yes",
        "-y/-n",
        help="Auto-accept downloads (default: yes). Use --no-yes to prompt.",
    ),
) -> None:
    """Install AiZynthFinder retrosynthesis tooling into modules/."""
    if yes:
        os.environ["HEDGEHOG_AUTO_INSTALL"] = "1"
    else:
        os.environ.pop("HEDGEHOG_AUTO_INSTALL", None)

    from hedgehog.setup import ensure_aizynthfinder

    project_root = Path(__file__).resolve().parents[2]
    config_path = ensure_aizynthfinder(project_root)
    console.print(f"[bold]AiZynthFinder installed.[/bold] Config: {config_path}")


@setup_app.command("nvmolkit-worker")
def setup_nvmolkit_worker(
    yes: bool = typer.Option(
        True,
        "--yes/--no-yes",
        "-y",
        help="Auto-accept downloads (no prompt).",
    ),
    python_bin: str | None = typer.Option(
        None,
        "--python",
        help="Python interpreter for worker venv (default: python3.12 -> 3.11 -> 3.10).",
    ),
) -> None:
    """Install isolated nvMolKit worker environment in .venv-nvmolkit-worker."""
    if yes:
        os.environ["HEDGEHOG_AUTO_INSTALL"] = "1"

    from hedgehog.setup import ensure_nvmolkit_worker

    project_root = Path(__file__).resolve().parents[2]
    worker_path = ensure_nvmolkit_worker(project_root, python_bin=python_bin)
    console.print(f"[bold]nvMolKit worker installed.[/bold] Entry: {worker_path}")


@setup_app.command("shepherd-worker")
def setup_shepherd_worker(
    yes: bool = typer.Option(
        False,
        "--yes",
        "-y",
        help="Auto-accept downloads (no prompt).",
    ),
    python_bin: str | None = typer.Option(
        None,
        "--python",
        help="Python interpreter for worker venv (default: python3.12 -> 3.11 -> 3.10).",
    ),
) -> None:
    """Install isolated Shepherd-Score worker environment in .venv-shepherd-worker."""
    if yes:
        os.environ["HEDGEHOG_AUTO_INSTALL"] = "1"

    from hedgehog.setup import ensure_shepherd_worker

    project_root = Path(__file__).resolve().parents[2]
    worker_path = ensure_shepherd_worker(project_root, python_bin=python_bin)
    console.print(f"[bold]Shepherd worker installed.[/bold] Entry: {worker_path}")


@app.command()
def tui(
    session: str | None = typer.Option(
        None,
        "--session",
        "-s",
        help="Resume TUI session by job id.",
    ),
) -> None:
    """
    Launch the interactive TUI (Text User Interface).

    Requires Node.js to be installed. The TUI provides a visual interface
    for configuring and running the pipeline.

    Examples
    --------
    \b
    uv run hedgehog tui
    uv run hedgehog tui -s job_123
    """
    import shutil
    import subprocess
    from pathlib import Path

    def _find_tui_dir() -> Path | None:
        """Find the `tui/` directory in a source checkout."""
        search_roots: list[Path] = []
        try:
            search_roots.append(Path.cwd().resolve())
        except OSError:
            pass
        try:
            search_roots.append(Path(__file__).resolve().parent)
        except OSError:
            pass

        for root in search_roots:
            for parent in [root] + list(root.parents):
                if (parent / "pyproject.toml").exists():
                    candidate = parent / "tui"
                    if candidate.exists():
                        return candidate
        return None

    tui_dir = _find_tui_dir()

    if tui_dir is None or not tui_dir.exists():
        console.print(
            "[red]Error:[/red] TUI directory not found. "
            "Please run from a source checkout that contains the `tui/` folder."
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
            raise typer.Exit(code=1) from e

    tui_env = os.environ.copy()
    if session:
        tui_env["HEDGEHOG_TUI_SESSION"] = session
    else:
        tui_env.pop("HEDGEHOG_TUI_SESSION", None)

    # Launch the TUI
    try:
        subprocess.run(
            ["node", str(dist_dir / "index.js")],
            cwd=tui_dir,
            check=True,
            env=tui_env,
        )
    except subprocess.CalledProcessError as e:
        console.print(f"[red]TUI exited with error:[/red] {e.returncode}")
        raise typer.Exit(code=e.returncode) from e
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    app()
