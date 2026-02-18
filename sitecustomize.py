from __future__ import annotations

import logging
import os

if os.environ.get("HEDGEHOG_SHOW_PANDASTOOLS_WARNINGS", "").strip() != "1":
    for logger_name in ("rdkit.Chem.PandasTools", "rdkit.Chem.PandasPatcher"):
        rdkit_logger = logging.getLogger(logger_name)
        rdkit_logger.setLevel(logging.ERROR)
        rdkit_logger.propagate = False
        if not any(
            isinstance(handler, logging.NullHandler)
            for handler in rdkit_logger.handlers
        ):
            rdkit_logger.addHandler(logging.NullHandler())
