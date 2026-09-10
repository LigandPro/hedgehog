"""Wave registry for structural filters stage.

Provides wave-based execution patterns for the struct_filters pipeline.
"""

# Wave definitions: ordered groups of filters that can run together.
# Each wave is a list of filter names; waves execute sequentially,
# filters within a wave may execute in parallel in the future.
DEFAULT_WAVES = [
    # Wave 1: fast filters (pure RDKit / SMARTS-based)
    [
        "common_alerts",
        "bredt",
        "protecting_groups",
        "ring_infraction",
        "stereo_center",
        "halogenicity",
        "symmetry",
    ],
    # Wave 2: heavier filters (medchem / external)
    ["molgraph_stats", "molcomplexity", "NIBR", "lilly"],
]


POLICY_CALCULATION_FILTER = {
    "undefined_stereo_center": "stereo_center",
}


def policy_calculation_filter(policy_name: str) -> str:
    """Return the calculated filter that supplies a policy pass mask."""
    return POLICY_CALCULATION_FILTER.get(policy_name, policy_name)


def get_calculated_policy_names(
    calculated_filters, config_struct_filters: dict
) -> list[str]:
    """Return report/enforcement policy names backed by calculated filters."""
    names = list(calculated_filters)
    if (
        "stereo_center" in names
        and "filter_undefined_stereo_center" in config_struct_filters
    ):
        index = names.index("stereo_center") + 1
        names.insert(index, "undefined_stereo_center")
    return names


def get_aligned_enforced_filters(
    config: dict, config_struct_filters: dict | None = None
) -> set[str] | None:
    """Resolve calculated filters that participate in Stage 3 survival.

    Per-filter ``filter_<name>`` booleans are the preferred policy. The legacy
    ``enforced_filters`` list remains supported when no filter flags are
    present. Alignment never overrides structural hard-gate policy.
    """
    del config  # kept for call-site compatibility; structural policy is config-only
    if config_struct_filters is not None:
        catalog_filters = [name for wave in DEFAULT_WAVES for name in wave]
        known_filters = [*catalog_filters, *POLICY_CALCULATION_FILTER]
        explicit_flags = {
            name: bool(config_struct_filters[f"filter_{name}"])
            for name in known_filters
            if f"filter_{name}" in config_struct_filters
        }
        if explicit_flags:
            calculated = {
                name
                for name in catalog_filters
                if config_struct_filters.get(f"calculate_{name}") is True
            }
            fully_explicit = calculated and all(
                f"filter_{name}" in config_struct_filters for name in calculated
            )
            if not calculated or fully_explicit:
                return {name for name, enabled in explicit_flags.items() if enabled}

            # Older configs calculated every listed method as a hard filter.
            # Preserve that behavior for missing flags while applying new flags
            # explicitly, so adding one policy toggle cannot disable unrelated
            # filters. The undefined-stereo toggle implicitly makes total stereo
            # diagnostic-only unless its own hard flag is explicitly enabled.
            enforced = set(calculated)
            for name, enabled in explicit_flags.items():
                if enabled:
                    enforced.add(name)
                else:
                    enforced.discard(name)
            if (
                "undefined_stereo_center" in explicit_flags
                and "stereo_center" not in explicit_flags
            ):
                enforced.discard("stereo_center")
            return enforced
        if "enforced_filters" in config_struct_filters:
            return _normalize_enforced_filters(
                config_struct_filters.get("enforced_filters") or []
            )
    return None


def _normalize_enforced_filters(raw_rules) -> set[str]:
    """Normalize rule-level entries to structural filter names."""
    enabled: set[str] = set()
    for raw_rule in raw_rules:
        rule = str(raw_rule)
        enabled.add("common_alerts" if rule.startswith("common_alerts:") else rule)
    return enabled
