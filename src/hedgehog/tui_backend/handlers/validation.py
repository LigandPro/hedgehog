"""Validation handler for TUI backend."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

from ..validators import ConfigValidator

if TYPE_CHECKING:
    from ..server import JsonRpcServer


class ValidationHandler:
    """Handler for validation-related RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server

    def validate_input_file(self, path: str) -> dict[str, Any]:
        """Validate that an input file exists and is readable."""
        file_path = Path(path).expanduser().resolve()

        result = {
            "valid": False,
            "path": str(file_path),
            "errors": [],
            "warnings": [],
        }

        if not file_path.exists():
            result["errors"].append(f"File not found: {path}")
            return result

        if not file_path.is_file():
            result["errors"].append(f"Not a file: {path}")
            return result

        if not file_path.suffix.lower() in [".sdf", ".csv", ".smi", ".smiles"]:
            result["warnings"].append(f"Unexpected file extension: {file_path.suffix}")

        # Check file size
        size = file_path.stat().st_size
        if size == 0:
            result["errors"].append("File is empty")
            return result

        if size > 1024 * 1024 * 1024:  # 1GB
            result["warnings"].append(
                "File is very large (>1GB), processing may be slow"
            )

        result["valid"] = len(result["errors"]) == 0
        result["size"] = size
        return result

    def validate_receptor_pdb(self, path: str) -> dict[str, Any]:
        """Validate a receptor PDB file for docking."""
        file_path = Path(path).expanduser().resolve()

        result = {
            "valid": False,
            "path": str(file_path),
            "errors": [],
            "warnings": [],
        }

        if not file_path.exists():
            result["errors"].append(f"Receptor file not found: {path}")
            return result

        if not file_path.is_file():
            result["errors"].append(f"Not a file: {path}")
            return result

        if file_path.suffix.lower() != ".pdb":
            result["warnings"].append(
                f"Expected .pdb extension, got: {file_path.suffix}"
            )

        # Basic PDB validation
        try:
            content = file_path.read_text()
            lines = content.split("\n")

            has_atom = any(
                line.startswith("ATOM") or line.startswith("HETATM") for line in lines
            )
            if not has_atom:
                result["errors"].append("No ATOM or HETATM records found in PDB file")
                return result

            atom_count = sum(1 for line in lines if line.startswith("ATOM"))
            hetatm_count = sum(1 for line in lines if line.startswith("HETATM"))

            result["atom_count"] = atom_count
            result["hetatm_count"] = hetatm_count

            if atom_count == 0:
                result["warnings"].append(
                    "No protein atoms found (only HETATM records)"
                )

        except Exception as e:
            result["errors"].append(f"Error reading PDB file: {e}")
            return result

        result["valid"] = len(result["errors"]) == 0
        return result

    def validate_output_directory(self, path: str) -> dict[str, Any]:
        """Validate output directory."""
        dir_path = Path(path).expanduser().resolve()

        result = {
            "valid": False,
            "path": str(dir_path),
            "errors": [],
            "warnings": [],
        }

        if dir_path.exists():
            if not dir_path.is_dir():
                result["errors"].append(f"Path exists but is not a directory: {path}")
                return result

            # Check if writable
            test_file = dir_path / ".hedgehog_test"
            try:
                test_file.write_text("test")
                test_file.unlink()
            except PermissionError:
                result["errors"].append(f"Directory is not writable: {path}")
                return result
            except Exception as e:
                result["errors"].append(f"Cannot write to directory: {e}")
                return result

            # Check if not empty
            if any(dir_path.iterdir()):
                result["warnings"].append(
                    "Output directory is not empty, files may be overwritten"
                )
        else:
            # Check if parent exists and is writable
            parent = dir_path.parent
            if not parent.exists():
                result["errors"].append(f"Parent directory does not exist: {parent}")
                return result

            result["warnings"].append("Directory does not exist and will be created")

        result["valid"] = len(result["errors"]) == 0
        return result

    def validate_config(
        self, config_type: str, config: dict[str, Any]
    ) -> dict[str, Any]:
        """Validate a configuration dictionary using centralized ConfigValidator."""
        return ConfigValidator.validate(config_type, config)
