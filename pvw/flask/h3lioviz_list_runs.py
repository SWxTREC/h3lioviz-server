"""Helpers for listing available H3lioviz runs."""

from __future__ import annotations

from decimal import Decimal
from typing import Any

import boto3

from helpers import _get_env_var

def list_runs() -> Any:
    table_name = _get_env_var("TABLE_NAME")

    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table(table_name)

    # Get all items from table using pagination
    data = []
    scan_kwargs = {}
    while True:
        response = table.scan(**scan_kwargs)
        data.extend(response.get("Items", []))
        if "LastEvaluatedKey" not in response:
            break
        scan_kwargs["ExclusiveStartKey"] = response["LastEvaluatedKey"]

    # Convert all Decimal types in table to floats for json conversion
    for item in data:
        for key in item:
            if isinstance(item[key], Decimal):
                item[key] = float(item[key])

    # Return the contents from metadata.json for all runs
    return data
