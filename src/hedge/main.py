import subprocess
from pathlib import Path
from typing import Optional
from enum import Enum

import matplotlib
matplotlib.use('Agg')

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')

from hedge.pipeline import calculate_metrics
from hedge.configs.logger import logger, load_config
from hedge.utils.data_prep import prepare_input_data

DEFAULT_CONFIG_PATH = './src/hedge/configs/config.yml'
SAMPLED_MOLS_FILENAME = 'sampledMols.csv'
STAGE_OVERRIDE_KEY = '_run_single_stage_override'


def _validate_input_path(input_path):
    """Validate input path and return Path object if valid, None otherwise."""
    if '*' in input_path or '?' in input_path:
        return None
    input_path_obj = Path(input_path)
    return input_path_obj if input_path_obj.exists() else None


def _get_unique_results_folder(base_folder):
    """
    Generate a unique folder name by appending a number if the folder already exists.

    If the base folder exists and contains any files or directories, creates a new folder
    with an incremented suffix (e.g., results -> results1 -> results2).

    Args:
        base_folder: Path object representing the base folder

    Returns:
        Path object with a unique folder name
    """
    base_folder = Path(base_folder)

    # If folder doesn't exist, return as is
    if not base_folder.exists():
        return base_folder

    # Check if folder exists and is not empty
    if base_folder.exists() and any(base_folder.iterdir()):
        # Extract base name and try to find a counter suffix
        base_name = base_folder.name
        parent = base_folder.parent

        # Try to extract existing counter from base name
        import re
        match = re.match(r'^(.+?)(\d+)$', base_name)
        if match:
            name_without_counter = match.group(1)
            start_counter = int(match.group(2)) + 1
        else:
            name_without_counter = base_name
            start_counter = 1

        # Find the next available folder name
        counter = start_counter
        while True:
            new_folder = parent / f"{name_without_counter}{counter}"
            if not new_folder.exists() or not any(new_folder.iterdir()):
                return new_folder
            counter += 1

    return base_folder


def preprocess_input_with_rdkit(input_path, folder_to_save, logger):
    """
    Preprocess input CSV file using RDKit (fallback when ligand_preparation_tool is not provided).
    Canonicalizes SMILES and removes duplicates within each model.
    
    Args:
        input_path: Path to the input CSV file
        folder_to_save: Output folder for prepared file
        logger: Logger instance
        
    Returns:
        Path to prepared CSV file if preprocessing was done, None otherwise
    """
    input_path_obj = _validate_input_path(input_path)
    if not input_path_obj:
        return None
    
    try:
        import pandas as pd
        from rdkit import Chem
        from rdkit.Chem import SanitizeMol
        
        folder_to_save.mkdir(parents=True, exist_ok=True)
        prepared_output = folder_to_save / f'prepared_{input_path_obj.stem}.csv'
        
        df = pd.read_csv(input_path_obj)
        cleaned_data = []
        
        for _, row in df.iterrows():
            smi = str(row['smiles'])
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                SanitizeMol(mol)
                canonical_smi = Chem.MolToSmiles(mol)
                row_data = {'smiles': canonical_smi}
                if 'model_name' in row:
                    row_data['model_name'] = row['model_name']
                cleaned_data.append(row_data)
            except Exception:
                continue
        
        if not cleaned_data:
            return None
        
        output_df = pd.DataFrame(cleaned_data)
        
        initial_count = len(output_df)
        output_df = output_df.drop_duplicates(subset=['smiles', 'model_name'], keep='first').reset_index(drop=True)
        duplicates_removed = initial_count - len(output_df)
        if duplicates_removed > 0:
            logger.info(f"Removed {duplicates_removed} duplicate molecules within models")
        
        output_df.to_csv(prepared_output, index=False)
        logger.info(f"RDKit preprocessing: {len(output_df)} molecules saved to {prepared_output}")
        return str(prepared_output)
    except Exception as e:
        logger.debug(f"RDKit preprocessing failed: {e}")
        return None


def preprocess_input_with_tool(input_path, ligand_preparation_tool, folder_to_save, logger):
    """
    Preprocess input file using external ligand preparation tool.
    
    Args:
        input_path: Path to the input file (CSV or SMI)
        ligand_preparation_tool: Path to the ligand preparation tool
        folder_to_save: Output folder for prepared file
        logger: Logger instance
        
    Returns:
        Path to prepared file if preprocessing was done, None otherwise
    """
    input_path_obj = _validate_input_path(input_path)
    if not input_path_obj:
        return None
    
    input_ext = input_path_obj.suffix.lower().lstrip('.')
    if input_ext == 'csv':
        input_format = '-icsv'
    elif input_ext in ('smi', 'ismi', 'cmi', 'txt'):
        input_format = '-ismi'
    else:
        return None
    
    folder_to_save.mkdir(parents=True, exist_ok=True)
    prepared_output = folder_to_save / f'prepared_{input_path_obj.stem}.csv'
    cmd = [ligand_preparation_tool, input_format, str(input_path), '-ucsv', str(prepared_output)]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return str(prepared_output) if prepared_output.exists() else None
    except Exception:
        return None

# Initialize Typer app and Rich console
app = typer.Typer(
    name="HEDGE",
    help="🦔 Hierarchical Evaluation of Drug GEnerators - Benchmark pipeline for generative models",
    add_completion=False,
    rich_markup_mode="rich"
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
    generated_mols_path: Optional[str] = typer.Option(None,
                                                      "--mols", "-m",
                                                      help="Path or glob pattern to generated SMILES files (overrides config)"
                                                     ),
    stage: Optional[Stage] = typer.Option(None,
                                          "--stage", "-s",
                                          help="Run only a specific pipeline stage. Please provide also --mols argument to specify the molecules path. If no --mols provided, the pipeline will use the molecules from the config file.",
                                          case_sensitive=False
                                         ),
    reuse_folder: bool = typer.Option(False,
                                      "--reuse",
                                      help="Force reuse of existing results folder (useful for reruns)"
                                     ),
    force_new_folder: bool = typer.Option(False,
                                          "--force-new",
                                          help="Force creation of a new results folder even when rerunning stages"
                                         ),
):
    """
    Run the molecular analysis pipeline.

    Examples:

    \b
    # Run full pipeline (auto-creates new folder if results exist)
    uv run hedge run

    \b
    # Rerun specific stage (auto-reuses existing folder)
    uv run hedge run --stage docking

    \b
    # Run stage with new molecules (auto-creates new folder)
    uv run hedge run --stage descriptors --mols data/*.csv

    \b
    # Force reuse existing folder
    uv run hedge run --reuse

    \b
    # Force create new folder even for stage rerun
    uv run hedge run --stage docking --force-new
    """

    # Display banner
    console.print(Panel.fit("[bold #B29EEE]🦔 HEDGE[/bold #B29EEE]\n"
                            "[dim]Hierarchical Evaluation of Drug GEnerators[/dim]\n"
                            "[dim italic]Developed by [bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim italic]",
                            border_style="#B29EEE"
                  ))
    
    config_dict = load_config(DEFAULT_CONFIG_PATH)
    
    if generated_mols_path:
        config_dict['generated_mols_path'] = generated_mols_path
        logger.info(f"[#B29EEE]Override:[/#B29EEE] Using molecules from: {generated_mols_path}")
    elif stage:
        logger.info(f"Using molecules from config: {config_dict['generated_mols_path']}")
    
    if stage:
        config_dict[STAGE_OVERRIDE_KEY] = stage.value
        logger.info(f"[#B29EEE]Override:[/#B29EEE] Running only stage: [bold]{stage.value}[/bold]")

    # Validate conflicting flags
    if reuse_folder and force_new_folder:
        logger.error("[red]Error:[/red] Cannot use --reuse and --force-new together. Please choose one.")
        raise typer.Exit(code=1)

    # Determine folder strategy based on flags and context
    original_folder = Path(config_dict['folder_to_save'])

    if reuse_folder:
        # Explicit reuse: always use configured folder
        folder_to_save = original_folder
        logger.info(f"[#B29EEE]Folder mode:[/#B29EEE] Reusing folder '{folder_to_save}' (--reuse flag)")
    elif force_new_folder:
        # Explicit new: always create incremented folder
        folder_to_save = _get_unique_results_folder(original_folder)
        if folder_to_save != original_folder:
            logger.info(f"[#B29EEE]Folder mode:[/#B29EEE] Creating new folder '{folder_to_save}' (--force-new flag)")
    else:
        # Automatic logic (Hybrid Variant 4)
        if stage and not generated_mols_path:
            # Stage rerun on existing data → reuse folder
            folder_to_save = original_folder
            logger.info(f"[#B29EEE]Folder mode:[/#B29EEE] Reusing folder '{folder_to_save}' for stage execution")
        else:
            # Full run OR stage with new molecules → create new folder
            folder_to_save = _get_unique_results_folder(original_folder)
            if folder_to_save != original_folder:
                logger.info(f"[#B29EEE]Folder mode:[/#B29EEE] Folder '{original_folder}' contains results. Using '{folder_to_save}' instead.")

    # Update config with the chosen folder path
    config_dict['folder_to_save'] = str(folder_to_save)

    ligand_preparation_tool = config_dict.get('ligand_preparation_tool')
    original_input_path = config_dict.get('generated_mols_path') or generated_mols_path
    if original_input_path:
        if ligand_preparation_tool:
            prepared_path = preprocess_input_with_tool(original_input_path, ligand_preparation_tool, folder_to_save, logger)
        else:
            prepared_path = preprocess_input_with_rdkit(original_input_path, folder_to_save, logger)
        
        if prepared_path:
            config_dict['generated_mols_path'] = prepared_path
            logger.info(f"Using preprocessed input: {prepared_path}")
    
    data = prepare_input_data(config_dict, logger)
    
    if 'mol_idx' not in data.columns or data['mol_idx'].isna().all():
        try:
            from hedge.utils.mol_index import assign_mol_idx
            data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)
        except Exception as e:
            logger.error(f'Failed to assign mol_idx: {e}')
            raise
    
    save_mols = config_dict.get('save_sampled_mols', False) or stage is not None
    if save_mols:
        folder_to_save.mkdir(parents=True, exist_ok=True)
        output_path = folder_to_save / SAMPLED_MOLS_FILENAME
        data.to_csv(output_path, index=False)
        logger.info(f"[#B29EEE]✓[/#B29EEE] Sampled total of {len(data)} molecules saved to {output_path}")
    
    # Run metrics calculation pipeline
    logger.info("[bold #B29EEE]Starting pipeline...[/bold #B29EEE]")
    calculate_metrics(data, config_dict)
    logger.info("[bold][#B29EEE]Ligand Pro[/#B29EEE] thanks you for using [#B29EEE]🦔 HEDGE[/#B29EEE]![/bold]")


@app.command()
def info():
    """
    Display pipeline information and available stages.
    """
    table = Table(title="Available Pipeline Stages", show_header=True, header_style="bold #B29EEE")
    table.add_column("Stage", style="#B29EEE", no_wrap=True)
    table.add_column("Description", style="white")
    
    stages_info = {"descriptors": "Compute 22 physicochemical descriptors per molecule",
                   "struct_filters": "Apply structural filters (Lilly, NIBR, PAINS, etc.)",
                   "synthesis": "Evaluate synthetic accessibility using retrosynthesis (AiZynthFinder) and other metrics",
                   "docking": "Calculate docking scores with Smina/Gnina"
                  }
    
    for stage_name, description in stages_info.items():
        table.add_row(stage_name, description)
    
    console.print(table)
    console.print("\n[dim]Example (1): uv run hedge run --stage descriptors[/dim]")
    console.print("[dim]Example (2): uv run hedge run --help [/dim]")



@app.command()
def version():
    """
    Display version information.
    """
    console.print("[bold #B29EEE]🦔 HEDGE[/bold #B29EEE] version [bold]1.0.0[/bold]")
    console.print("[dim]Hierarchical Evaluation of Drug GEnerators[/dim]")
    console.print("[dim]Developed by [bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim]")


if __name__ == '__main__':
    app()
    