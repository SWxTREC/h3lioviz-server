# Flask

## Installing Additional Dependencies

Use [../pvw/requirements.txt](../pvw/requirements.txt) to add any needed dependencies you want Docker to install.

## Flask API Endpoints

The Flask API Endpoints are supposed to mirror those in the swxtrec-cdk, but may lag behind as the Flask routes are not in active use.
Instead we are hitting h3lioviz-api.\* which is an API Gateway that routes to AWS Lambda Functions.

The Flask server provides the following REST endpoints at localhost:5000:

- `GET /h3lioviz/metadata/health` - Health check endpoint, returns HTTP 200 code with body: `{"status": "ok"}`
- `GET /h3lioviz/metadata/getTimeSeries/<run_id>/<satellite>` - Retrieves time-series data for a specific run and satellite
- `GET /h3lioviz/metadata/availableRuns` - Lists all available simulation runs
- `GET /h3lioviz/metadata/syncMetadata` - Synchronizes metadata from data sources

## Logging

The logs include the following:

- Server startup information
- Any exceptions that occurred
- Basic information on all requests including IP address, request method, date, browser, and status code

  `/data/launcher/log/flask.log`
