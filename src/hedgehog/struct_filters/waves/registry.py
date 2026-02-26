"""Wave registry for structural filters stage.

Provides wave-based execution patterns for the struct_filters pipeline.
"""

from hedgehog.configs.logger import logger

# Wave definitions: ordered groups of filters that can run together.
# Each wave is a list of filter names; waves execute sequentially,
# filters within a wave may execute in parallel in the future.
DEFAULT_WAVES = [
    # Wave 1: fast filters (pure RDKit / SMARTS-based)
    ["common_alerts", "bredt", "protecting_groups", "ring_infraction",
     "stereo_center", "halogenicity", "symmetry"],
    # Wave 2: heavier filters (medchem / external)
    ["molgraph_stats", "molcomplexity", "NIBR", "lilly"],
]


def get_waves(config_struct_filters=None):
    """Return wave definitions, optionally filtered by enabled filters.

    Args:
        config_struct_filters: If provided, only include filters that are
            enabled (``calculate_<name>: true``) in the config.

    Returns:
        list[list[str]]: Waves of filter names.
    """
    waves = DEFAULT_WAVES

    if config_struct_filters is not None:
        enabled = {
            k.replace("calculate_", "")
            for k, v in config_struct_filters.items()
            if "calculate_" in k and v
        }
        waves = [
            [f for f in wave if f in enabled]
            for wave in waves
        ]
        waves = [w for w in waves if w]

    return waves
