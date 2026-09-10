"""Small shared helpers for structural filters config normalization."""


def is_include_rulesets_all(raw) -> bool:
    """Return True when include_rulesets means every catalog ruleset."""
    return isinstance(raw, str) and raw.strip().lower() == "all"


def resolve_include_rulesets(raw, available_rulesets) -> list[str]:
    """Normalize ``include_rulesets`` to a concrete list of ruleset names.

    Accepted forms:
    - ``None`` / missing / ``[]`` → no rulesets
    - ``"all"`` → every name in *available_rulesets* (catalog order preserved)
    - ``[name, ...]`` → those names

    Raises:
        ValueError: on unsupported types or an explicit ``all`` entry inside a list.
    """
    if raw is None:
        return []

    if is_include_rulesets_all(raw):
        return [str(name) for name in available_rulesets]

    if isinstance(raw, str):
        raise ValueError(
            "include_rulesets must be 'all', null/[], or a list of ruleset names; "
            f"got string {raw!r}."
        )

    if not isinstance(raw, (list, tuple, set)):
        raise ValueError(
            "include_rulesets must be 'all', null/[], or a list of ruleset names; "
            f"got {type(raw).__name__}."
        )

    names: list[str] = []
    for value in raw:
        if value is None:
            continue
        name = str(value).strip()
        if not name:
            continue
        if name.lower() == "all":
            raise ValueError(
                "include_rulesets: use scalar 'all' for every ruleset, "
                "not a list entry named 'all'."
            )
        names.append(name)
    return names
