import os
import polars as pl
from pathlib import Path

from hedge.configs.logger import logger, load_config
from hedge.stages.synthesis.utils import *

def _order_identity_columns(df):
    """Order dataframe columns with identity columns first."""
    id_cols = ['smiles', 'model_name', 'mol_idx']
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered_cols = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df.select(ordered_cols)


def main(config):
    """Main entry point for synthesis stage (retrosynthesis analysis).

    Args:
        config: Configuration dictionary containing pipeline settings
    """
    folder_to_save = process_path(config['folder_to_save'])
    output_folder_str = os.path.join(folder_to_save, 'stages', '04_synthesis')
    output_folder = Path(output_folder_str)
    output_folder.mkdir(parents=True, exist_ok=True)
    config_synthesis = load_config(config['config_synthesis'])
    input_path = get_input_path(config, folder_to_save)
    
    if not os.path.exists(input_path):
        logger.error(f"Input file not found: {input_path}")
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    try:
        input_df = pl.read_csv(input_path)
        logger.info(f"Processing {len(input_df)} filtered molecules")
    except Exception as e:
        logger.error(f"Could not load input data: {e}")
        raise
    
    if len(input_df) == 0:
        logger.warning("No molecules to process for synthesis analysis")
        return
    
    scored_df = calculate_synthesis_scores(input_df, folder_to_save, config_synthesis)
    scores_output = output_folder / 'synthesis_scores.csv'
    scored_df_ordered = _order_identity_columns(scored_df)
    scored_df_ordered.write_csv(scores_output)
    score_filtered_df = apply_synthesis_score_filters(scored_df, config_synthesis)

    if len(score_filtered_df) == 0:
        logger.warning("No molecules passed synthesis score filters")
        output_file = output_folder / 'filtered_molecules.csv'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        score_filtered_df.write_csv(output_file)
        logger.info(f"Saved 0 molecules to {output_file}")
        return
    
    input_smiles_file = output_folder / 'input_smiles.smi'
    retrosynth_module = Path(__file__).parent.parent.parent.parent.parent / 'modules' / 'retrosynthesis' / 'aizynthfinder'
    aizynth_config_file = retrosynth_module / 'public' / 'config.yml'
    
    if not aizynth_config_file.exists():
        logger.error(f"AiZynthFinder config file not found: {aizynth_config_file}")
        logger.error(f"To set up retrosynthesis, run:")
        logger.error(f"  cd {retrosynth_module}")
        logger.error(f"  mkdir -p public")
        logger.error(f"  uv run python -m aizynthfinder.tools.download_public_data ./public")
        logger.error(f"Synthesis stage will be skipped - continuing pipeline without retrosynthesis")
        output_file = output_folder / 'filtered_molecules.csv'
        score_filtered_df_copy = _order_identity_columns(score_filtered_df)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        score_filtered_df_copy.write_csv(output_file)
        return

    prepare_input_smiles(score_filtered_df, input_smiles_file)
    output_json_file = output_folder / 'retrosynthesis_results.json'
    success = run_aizynthfinder(input_smiles_file=input_smiles_file,
                                output_json_file=output_json_file,
                                config_file=aizynth_config_file
                            )
    if not success:
        logger.error("Retrosynthesis analysis failed")
        raise RuntimeError("Retrosynthesis analysis failed")
    
    retrosynth_df = parse_retrosynthesis_results(output_json_file)
    
    if len(retrosynth_df) == 0:
        logger.warning("No retrosynthesis results found in JSON file")
        return
    
    merged_df = merge_retrosynthesis_results(score_filtered_df, retrosynth_df)
    extended_output = output_folder / 'synthesis_extended.csv'
    merged_df_ordered = _order_identity_columns(merged_df)
    merged_df_ordered.write_csv(extended_output)
    filter_solved_only = config_synthesis.get('filter_solved_only', True)

    if filter_solved_only:
        filtered_df = merged_df.filter(pl.col('solved') == 1)
    else:
        filtered_df = merged_df.clone()
        logger.info("Keeping all molecules (filter_solved_only=False)")

    output_file = output_folder / 'filtered_molecules.csv'
    filtered_df = _order_identity_columns(filtered_df)
    filtered_df.write_csv(output_file)
    
    avg_time = merged_df['search_time'].mean() if 'search_time' in merged_df.columns else 0.0
    logger.info(f"Average retrosynthesis search time: {avg_time:.2f} s")