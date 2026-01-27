"""File handling for TUI backend."""
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..server import JsonRpcServer


class FilesHandler:
    """Handler for file-related RPC methods."""

    def __init__(self, server: 'JsonRpcServer'):
        self.server = server

    def list_files(self, path: str, extensions: list[str] | None = None) -> list[str]:
        """List files in a directory, optionally filtered by extension."""
        dir_path = Path(path).expanduser().resolve()

        if not dir_path.exists():
            raise FileNotFoundError(f'Directory not found: {path}')

        if not dir_path.is_dir():
            raise NotADirectoryError(f'Not a directory: {path}')

        # Normalize extensions: handle both '.csv' and 'csv' formats
        normalized_extensions = None
        if extensions is not None:
            normalized_extensions = {ext.lstrip('.').lower() for ext in extensions}

        files = []
        for item in sorted(dir_path.iterdir()):
            if item.is_file():
                if normalized_extensions is None:
                    files.append(str(item))
                elif item.suffix.lstrip('.').lower() in normalized_extensions:
                    files.append(str(item))

        return files

    def list_directory(self, path: str, extensions: list[str] | None = None) -> list[dict[str, Any]]:
        """List directory contents with metadata."""
        dir_path = Path(path).expanduser().resolve()

        if not dir_path.exists():
            raise FileNotFoundError(f'Directory not found: {path}')

        if not dir_path.is_dir():
            raise NotADirectoryError(f'Not a directory: {path}')

        # Normalize extensions: handle both '.csv' and 'csv' formats
        normalized_extensions = None
        if extensions is not None:
            normalized_extensions = {ext.lstrip('.').lower() for ext in extensions}

        entries = []

        try:
            for item in sorted(dir_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
                # Skip hidden files
                if item.name.startswith('.'):
                    continue

                is_dir = item.is_dir()

                # Filter by extension for files
                if not is_dir and normalized_extensions is not None:
                    file_ext = item.suffix.lstrip('.').lower()
                    if file_ext not in normalized_extensions:
                        continue
                
                entry = {
                    'name': item.name,
                    'path': str(item),
                    'isDirectory': is_dir,
                }
                
                if not is_dir:
                    try:
                        stat = item.stat()
                        entry['size'] = stat.st_size
                        entry['modified'] = stat.st_mtime
                    except OSError:
                        pass
                
                entries.append(entry)
        except PermissionError:
            raise PermissionError(f'Permission denied: {path}')
        
        return entries

    def count_molecules(self, path: str) -> dict[str, Any]:
        """Count molecules in a file (CSV, SMI, SDF)."""
        file_path = Path(path).expanduser().resolve()

        if not file_path.exists():
            raise FileNotFoundError(f'File not found: {path}')

        if not file_path.is_file():
            raise ValueError(f'Not a file: {path}')

        ext = file_path.suffix.lower()
        count = 0

        try:
            if ext in ('.csv',):
                # CSV: count lines minus header
                with open(file_path, 'r') as f:
                    count = sum(1 for _ in f) - 1  # subtract header
                    count = max(0, count)
            elif ext in ('.smi', '.smiles'):
                # SMI: count non-empty lines
                with open(file_path, 'r') as f:
                    count = sum(1 for line in f if line.strip())
            elif ext in ('.sdf', '.mol2'):
                # SDF: count $$$$ delimiters
                with open(file_path, 'r') as f:
                    content = f.read()
                    count = content.count('$$$$')
                    if count == 0:
                        # Fallback: might be single molecule
                        count = 1 if content.strip() else 0
            else:
                # Unknown format - try counting lines
                with open(file_path, 'r') as f:
                    count = sum(1 for line in f if line.strip())
        except Exception as e:
            return {'count': 0, 'error': str(e)}

        return {'count': count, 'path': str(file_path)}
