import os
from pathlib import Path
from typing import Optional
from enum import Enum

import matplotlib
import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')

matplotlib.use('Agg')

from calculate_metrics import calculate_metrics
from logger_config import logger
from configs.config_utils import load_config
from utils.data_prep import prepare_input_data

DEFAULT_CONFIG_PATH = './configs/config.yml'
SAMPLED_MOLS_FILENAME = 'sampledMols.csv'
STAGE_OVERRIDE_KEY = '_run_single_stage_override'

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
):
    """
    Run the molecular analysis pipeline.
    
    Examples:
    
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
    console.print(Panel.fit("[bold #B29EEE]🦔 HEDGE[/bold #B29EEE]\n"
                            "[dim]Hierarchical Evaluation of Drug GEnerators[/dim]\n"
                            "[dim italic]Developed by [bold #B29EEE]Ligand Pro[/bold #B29EEE][/dim italic]",
                            border_style="#B29EEE"
                  ))
    
    # Load configuration
    config_dict = load_config(DEFAULT_CONFIG_PATH)
    
    # Apply command-line overrides
    if generated_mols_path:
        config_dict['generated_mols_path'] = generated_mols_path
        logger.info(f"[yellow]Override:[/yellow] Using molecules from: {generated_mols_path}")
    # If stage is specified without --mols, use molecules from config
    elif stage:
        logger.info(f"Using molecules from config: {config_dict['generated_mols_path']}")
    
    if stage:
        config_dict[STAGE_OVERRIDE_KEY] = stage.value
        logger.info(f"[yellow]Override:[/yellow] Running only stage: [bold]{stage.value}[/bold]")
    
    # Prepare input data
    folder_to_save = Path(config_dict['folder_to_save'])
    data = prepare_input_data(config_dict, logger)
    
    # Always save sampled molecules when running a stage independently
    save_mols = config_dict.get('save_sampled_mols', False) or stage is not None
    if save_mols:
        folder_to_save.mkdir(parents=True, exist_ok=True)
        output_path = folder_to_save / SAMPLED_MOLS_FILENAME
        data.to_csv(output_path, index=False)
        logger.info(f"[green]✓[/green] Sampled {len(data)} molecules saved to {output_path}")
    
    # Run metrics calculation pipeline
    logger.info("[bold green]Starting pipeline...[/bold green]")
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
                   "docking": "Calculate docking scores with Smina/Gnina"
                  }
    
    for stage_name, description in stages_info.items():
        table.add_row(stage_name, description)
    
    console.print(table)
    console.print("\n[dim]Example: uv run python main.py run --stage descriptors[/dim]")


@app.command()
def version():
    """
    Display version information.
    """
    console.print("[bold #B29EEE]🦔 HEDGE[/bold #B29EEE] version [bold]1.0.0[/bold]")
    console.print("[dim]Hierarchical Evaluation of Drug GEnerators[/dim]")
    console.print("[dim]Developed by [bold]Ligand Pro[/bold][/dim]")


if __name__ == '__main__':
    app()