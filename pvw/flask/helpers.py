import os
from decimal import Decimal
from pathlib import Path
from typing import Any

DATA_PATH = Path("/data")


def _convert_decimals(obj: Any) -> Any:
    if isinstance(obj, list):
        return [_convert_decimals(item) for item in obj]
    if isinstance(obj, dict):
        return {key: _convert_decimals(value) for key, value in obj.items()}
    if isinstance(obj, Decimal):
        return float(obj)
    return obj


def _get_env_var(var_name: str) -> str:
    """
    Gets the value of the given environment variable. Empty variables are treated as unset and will raise an exception.
    """
    env_var = os.environ.get(var_name)
    if not env_var:
        raise RuntimeError(f"Missing required environment variable: {var_name}")

    return env_var
