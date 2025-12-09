# Flask

## Installing Additional Dependencies

Use [../pvw/requirements.txt](../pvw/requirements.txt) to add any needed dependencies you want Docker to install.

## Flask API Endpoints

The Flask server provides the following REST endpoints at localhost:5000:

- `GET /h3lioviz/metadata/health` - Health check endpoint, returns `{"status": "ok"}`
- `GET /h3lioviz/metadata/getTimeSeries/<run_id>/<satellite>` - Retrieves time-series data for a specific run and satellite
- `GET /h3lioviz/metadata/availableRuns` - Lists all available simulation runs
- `GET /h3lioviz/metadata/syncMetadata` - Synchronizes metadata from data sources

## Logging

`/data/launcher/log/flask.log`