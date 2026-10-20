# Flask

## Installing Additional Dependencies

Use [../pvw/requirements.txt](../pvw/requirements.txt) to add any needed dependencies you want Docker to install.

## Flask API Endpoints

The Flask API Endpoints are supposed to mirror those in the swxtrec-cdk, but may lag behind as the frontend does not call the Flask metadata routes.
Instead the frontend calls h3lioviz-api.\*, a CloudFront distribution that forwards to a load balancer, which routes to AWS Lambda Functions.
The exception is the health route: the load balancer uses it as its health check for the EC2.

The Flask server listens on 0.0.0.0:5000 inside the container, and Apache proxies `/h3lioviz/metadata/` to it. It provides the following REST endpoints:

- `GET /h3lioviz/metadata/health` - Health check endpoint, returns HTTP 200 code with body: `{"status": "ok"}`
- `GET /h3lioviz/metadata/getTimeSeries/<run_id>/<satellite>` - Retrieves time-series data for a specific run and satellite from the run's `evo.<satellite>` file in S3. `<satellite>` is the file name after `evo.`, for example `earth.json` (a trailing `d`, as in `earth.jsond`, is removed). Needs S3_BUCKET_NAME.
- `GET /h3lioviz/metadata/availableRuns` - Lists all available simulation runs from the DynamoDB table. Needs TABLE_NAME and AWS_DEFAULT_REGION.
- `GET /h3lioviz/metadata/syncMetadata` - Synchronizes the DynamoDB table with the runs in S3: adds runs it finds in S3 that have no table entry, and deletes table entries for runs it does not find in S3. Note that this GET request changes data. Needs S3_BUCKET_NAME, TABLE_NAME and AWS_DEFAULT_REGION.

A route whose environment variable is missing returns HTTP 500.

## Logging

The logs include the following:

- Server startup information
- Any exceptions that occurred
- Basic information on all requests including IP address, request method, date, and status code. The IP address is Apache's, since requests arrive through its proxy.

  `/data/launcher/log/flask.log`
