"""Helpers for listing available H3lioviz runs."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

DATA_PATH = Path("/data")


def list_runs() -> Any:

    data = []
    for path in DATA_PATH.iterdir():
        match = re.search(r"^pv-ready-data-(.*)$", path.name)
        if not match:
            continue

        with open(path / "metadata.json") as fp:
            data.append(json.load(fp))

    return data
