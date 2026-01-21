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
        
        files = []
        for item in sorted(dir_path.iterdir()):
            if item.is_file():
                if extensions is None or item.suffix.lstrip('.') in extensions:
                    files.append(str(item))
        
        return files

    def list_directory(self, path: str, extensions: list[str] | None = None) -> list[dict[str, Any]]:
        """List directory contents with metadata."""
        dir_path = Path(path).expanduser().resolve()
        
        if not dir_path.exists():
            raise FileNotFoundError(f'Directory not found: {path}')
        
        if not dir_path.is_dir():
            raise NotADirectoryError(f'Not a directory: {path}')
        
        entries = []
        
        try:
            for item in sorted(dir_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
                # Skip hidden files
                if item.name.startswith('.'):
                    continue
                    
                is_dir = item.is_dir()
                
                # Filter by extension for files
                if not is_dir and extensions is not None:
                    if item.suffix.lstrip('.') not in extensions:
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
