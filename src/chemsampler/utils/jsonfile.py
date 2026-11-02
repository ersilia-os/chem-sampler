"""Reading the small JSON settings files: a project's, a config file, `hpc.json`."""

import json
from pathlib import Path


def is_object(value) -> bool:
    return isinstance(value, dict)


def is_text(value) -> bool:
    return isinstance(value, str)


def is_list_of(value, kind) -> bool:
    return isinstance(value, list) and all(isinstance(item, kind) for item in value)


def read_object(path, *, expected: str = "") -> dict:
    """
    Read a JSON file that must hold an object.

    Parameters
    ----------
    path : str or Path
        The file to read.
    expected : str, optional
        Added to the error when the file holds something else, e.g. the keys it
        should have.

    Raises
    ------
    ValueError
        If the file is missing, is not valid JSON, or does not hold an object.
    """
    path = Path(path)
    try:
        text = path.read_text()
    except (FileNotFoundError, NotADirectoryError):
        raise ValueError(f"{path} does not exist") from None
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: not valid JSON ({e})") from None
    if not is_object(data):
        raise ValueError(f"{path}: expected a JSON object{expected}")
    return data
