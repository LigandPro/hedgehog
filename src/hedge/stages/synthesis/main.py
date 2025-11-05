from pathlib import Path
from typing import Any

import pandas as pd

from hedge.configs.logger import load_config, logger
from hedge.stages.synthesis.utils import (
    apply_synthesis_score_filters,
    calculate_synthesis_scores,
    get_input_path,
    merge_retrosynthesis_results,
    parse_retrosynthesis_results,
    prepare_input_smiles,
    process_path,
    run_aizynthfinder,
)


def _order_identity_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Order dataframe columns with identity columns first."""
    id_cols = ["smiles", "model_name", "mol_idx"]
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered_cols = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df[ordered_cols]


def main(config: dict[str, Any]) -> None:  # noqa: PLR0915
    """Run synthesis stage (retrosynthesis analysis).

    Args:
        config: Configuration dictionary containing pipeline settings
    """
    folder_to_save = process_path(config["folder_to_save"])
    output_folder_str = process_path(folder_to_save, "Synthesis")
    output_folder = Path(output_folder_str)
    config_synthesis = load_config(config["config_synthesis"])
    input_path = get_input_path(config, folder_to_save)

    if not Path(input_path).exists():
        msg = f"Input file not found: {input_path}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    try:
        input_df = pd.read_csv(input_path)
        logger.info(f"Processing {len(input_df)} filtered molecules")
    except Exception as e:
        logger.error(f"Could not load input data: {e}")
        raise

    if len(input_df) == 0:
        logger.warning("No molecules to process for synthesis analysis")
        return

    scored_df = calculate_synthesis_scores(input_df, folder_to_save, config_synthesis)
    scores_output = output_folder / "synthesis_scores.csv"
    scored_df_ordered = _order_identity_columns(scored_df)
    scored_df_ordered.to_csv(scores_output, index=False)
    score_filtered_df = apply_synthesis_score_filters(scored_df, config_synthesis)

    if len(score_filtered_df) == 0:
        logger.warning("No molecules passed synthesis score filters")
        output_file = output_folder / "passSynthesisSMILES.csv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        score_filtered_df.to_csv(output_file, index=False)
        logger.info(f"Saved 0 molecules to {output_file}")
        return

    input_smiles_file = output_folder / "input_smiles.smi"
    retrosynth_module = (
        Path(__file__).parent.parent.parent.parent.parent
        / "modules"
        / "retrosynthesis"
        / "aizynthfinder"
    )
    aizynth_config_file = retrosynth_module / "public" / "config.yml"

    if not aizynth_config_file.exists():
        logger.error(f"AiZynthFinder config file not found: {aizynth_config_file}")
        logger.error("To set up retrosynthesis, run:")
        logger.error(f"  cd {retrosynth_module}")
        logger.error("  mkdir -p public")
        logger.error(
            "  uv run python -m aizynthfinder.tools.download_public_data ./public"
        )
        logger.error(
            "Synthesis stage will be skipped - continuing pipeline without"
            " retrosynthesis"
        )
        output_file = output_folder / "passSynthesisSMILES.csv"
        score_filtered_df_copy = _order_identity_columns(score_filtered_df)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        score_filtered_df_copy.to_csv(output_file, index=False)
        return

    prepare_input_smiles(score_filtered_df, input_smiles_file)
    output_json_file = output_folder / "retrosynthesis_results.json"
    success = run_aizynthfinder(
        input_smiles_file=input_smiles_file,
        output_json_file=output_json_file,
        config_file=aizynth_config_file,
    )
    if not success:
        msg = "Retrosynthesis analysis failed"
        logger.error(msg)
        raise RuntimeError(msg)

    retrosynth_df = parse_retrosynthesis_results(output_json_file)

    if len(retrosynth_df) == 0:
        logger.warning("No retrosynthesis results found in JSON file")
        return

    merged_df = merge_retrosynthesis_results(score_filtered_df, retrosynth_df)
    extended_output = output_folder / "synthesis_extended.csv"
    merged_df_ordered = _order_identity_columns(merged_df)
    merged_df_ordered.to_csv(extended_output, index=False)
    filter_solved_only = config_synthesis.get("filter_solved_only", True)

    if filter_solved_only:
        filtered_df = merged_df[merged_df["solved"] == 1].copy()
    else:
        filtered_df = merged_df.copy()
        logger.info("Keeping all molecules (filter_solved_only=False)")

    output_file = output_folder / "passSynthesisSMILES.csv"
    filtered_df = _order_identity_columns(filtered_df)
    filtered_df.to_csv(output_file, index=False)

    if "search_time" in merged_df.columns:
        avg_time = merged_df["search_time"].mean()
    else:
        avg_time = 0.0
    logger.info(f"Average retrosynthesis search time: {avg_time:.2f} s")
