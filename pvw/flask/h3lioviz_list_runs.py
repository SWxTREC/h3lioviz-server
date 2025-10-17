"""Helpers for listing available H3lioviz runs."""

from __future__ import annotations

from decimal import Decimal
from typing import Any
import os

import boto3

def list_runs() -> Any:
    table_name = os.environ["TABLE_NAME"]

    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table(table_name)

    # Get all items from table
    data = table.scan()["Items"]

    # Convert all Decimal types in table to floats for json conversion
    for item in data:
        for key in item:
            if isinstance(item[key], Decimal):
                item[key] = float(item[key])

    # Return the contents from metadata.json for all runs
    return data
