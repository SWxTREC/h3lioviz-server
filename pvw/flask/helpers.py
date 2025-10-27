from typing import Any
from decimal import Decimal
import os


def _convert_decimals(obj: Any) -> Any:
    if isinstance(obj, list):
        return [_convert_decimals(item) for item in obj]
    if isinstance(obj, dict):
        return {key: _convert_decimals(value) for key, value in obj.items()}
    if isinstance(obj, Decimal):
        return float(obj)
    return obj

def _get_env_var(var_name: str) -> str:
    env_var = os.environ.get(var_name)
    if not env_var:
        raise RuntimeError(f"Missing required environment variable: {var_name}")

    return env_var
