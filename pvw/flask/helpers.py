import os
from typing import Any
from decimal import Decimal


def _require_env(var_name: str) -> str:
    value = os.environ.get(var_name)
    if not value:
        raise RuntimeError(f"Missing required environment variable: {var_name}")
    return value


def _convert_decimals(obj: Any) -> Any:
    if isinstance(obj, list):
        return [_convert_decimals(item) for item in obj]
    if isinstance(obj, dict):
        return {key: _convert_decimals(value) for key, value in obj.items()}
    if isinstance(obj, Decimal):
        return float(obj)
    return obj