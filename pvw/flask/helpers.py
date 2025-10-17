from typing import Any
from decimal import Decimal

def _convert_decimals(obj: Any) -> Any:
    if isinstance(obj, list):
        return [_convert_decimals(item) for item in obj]
    if isinstance(obj, dict):
        return {key: _convert_decimals(value) for key, value in obj.items()}
    if isinstance(obj, Decimal):
        return float(obj)
    return obj