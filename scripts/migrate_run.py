# Purpose: This script was developed because we no longer have the original data
#          for some runs and therefore need to apply updates to runs that have
#          already gone through some processing.
#
# Note: This script does not update metadata.json for the v0.0.0 -> v1.1.1 migration
#       besides updating/adding h3lioviz_processing_version. This is done because some
#       data in metadata.json for old runs appears to have been added manually (like the institute)
#
import argparse
import json
import logging
import shutil
from pathlib import Path
from datetime import datetime

import xarray as xr

logger = logging.getLogger(__name__)


def migrate_run_v0_to_v1(path: Path):
    # This is intended to fail if the path already exists to prevent overwriting
    output_path = Path(f"{path!s}-migrated")
    output_path.mkdir()

    tim_paths = path.glob("pv-tim.*.nc")
    for tim_path in tim_paths:
        with xr.load_dataset(tim_path) as ds:
            ds = migrate_tim_v0_to_v1(ds)

            ds.to_netcdf(
                output_path / tim_path.name,
                encoding={"time": {"units": "seconds since 1970-01-01"}},
            )

    evo_paths = path.glob("evo.*.nc")
    for evo_path in evo_paths:
        new_path = output_path / evo_path.name
        with xr.load_dataset(evo_path) as ds:
            ds = migrate_evo_v0_to_v1(ds)

            ds.to_netcdf(
                new_path,
                encoding={"time": {"units": "seconds since 1970-01-01"}},
            )

        with xr.load_dataset(new_path) as ds:
            json_file = Path(str(new_path).replace(".nc", ".json"))

            # Save the json file in the same way as process_output.py
            with open(json_file, "w") as fp:
                fp.write(
                    json.dumps(
                        json.loads(
                            json.dumps(ds.to_dict(), default=serialize_datetime),
                            parse_float=lambda x: round(float(x), 3),
                        )
                    )
                )

    # Add/update h3lioviz_processing_version in metadata.json
    with open(path / "metadata.json", "r") as fp:
        metadata = json.load(fp)
        metadata["h3lioviz_processing_version"] = "1.0.0"
    with open(output_path / "metadata.json", "w") as fp:
        json.dump(metadata, fp)

    # Copy the cone file unmodified if it exists
    cone_file_path = next(path.glob("cone2bc.in.*"), None)
    if cone_file_path:
        shutil.copyfile(cone_file_path, output_path / cone_file_path.name)


def migrate_tim_v0_to_v1(ds: xr.Dataset):
    # Updated as a part of "DB-3308" PR Oct. 24, 2025
    ds["Pressure"] = ds["Pressure"] / 1e6
    ds["Pressure"].attrs.update({"units": "r^2 * N / m^3 * km^2 / s^2"})

    # Updates as a part of "Multiply Br by normalized radius squared" PR on Oct. 21, 2025
    ds["Br"] *= ds["radius"] ** 2

    return ds


# Note that the Br = B1 is the only thing updated. However during normal processing B1 is
# dropped so the odds that any of the evo files are updated are low for this migration.
def migrate_evo_v0_to_v1(ds: xr.Dataset):
    # Updated as a part of "Added Br to evo" PR on Oct. 13, 2025
    if "B1" in ds:
        ds["Br"] = ds["B1"]
    else:
        logger.warning("Variable B1 does not exist in dataset.")

    return ds


def serialize_datetime(obj):
    if isinstance(obj, datetime):
        return obj.timestamp()
    else:
        return obj


def main():
    parser = argparse.ArgumentParser(prog="Migrate H3lioviz Run")

    parser.add_argument(
        "path",
        type=Path,
        help="Path to the directory containing .nc files that need migration.",
    )

    args = parser.parse_args()

    if not args.path.exists() or args.path.is_file():
        raise ValueError(
            f"Provided path '{args.path!s}' does not exist or is not a directory."
        )

    migrate_run_v0_to_v1(args.path)


if __name__ == "__main__":
    main()
