"""Custom exception hierarchy for HEDGEHOG pipeline."""


class HedgehogError(Exception):
    """Base exception for all HEDGEHOG errors."""


class StageError(HedgehogError):
    """Error during pipeline stage execution."""


class DockingError(StageError):
    """Error during molecular docking."""


class FilterError(StageError):
    """Error during filtering operations."""


class ConfigError(HedgehogError):
    """Error in configuration loading or validation."""


class MoleculeError(HedgehogError):
    """Error in molecule processing (RDKit, SMILES, etc.)."""
