# H3lioViz

This repository hosts the code to run a ParaView visualization server for 3D heliospheric output from codes such as Enlil and Euhforia, as well as a way to fetch the associated metadata and host the frontend.

Live deployment for interactive use: <https://swx-trec.com/h3lioviz/>

## Additional Documentation

- [Apache](docs/Apache.md)
- [Flask](docs/Flask.md)
- [Paraview](docs/Paraview.md)

## Overview

The H3lioViz server is a containerized application that provides:

- **ParaView Web Service**: 3D visualization of heliospheric simulation data
- **Flask Metadata API**: Endpoints for retrieving run information and time-series data
- **Apache Web Server**: Serves the frontend and acts as a reverse proxy
- **ParaView Pre-Processing**: Scripts that generate the data used by the ParaView Web Service (not containerized). Found in the `scripts/` directory.

### Architecture

Apache serves as the main entry point and handles routing:

- `/` -> Reroutes to `/h3lioviz/`
- `/h3lioviz/*` -> All routes except the below endpoints route to the Frontend Web Application
- `/h3lioviz/paraview` -> ParaView Web service (creates visualization sessions)
- `/h3lioviz/proxy` -> WebSocket proxy (maps session IDs to ParaView ports)
- `/h3lioviz/metadata/` -> Flask API for metadata operations

### Container Startup Script

The docker/scripts/server.sh script initializes the docker container by using certain environment variables to update configuration files.

server.sh (paraview/websockets) environment variables:

- SERVER_NAME: the server name to use for the session url. Gets returned from paraview-web to tell the frontend where to connect to the websocket. Note that the current routing expects {domain}/h3lioviz.
- PROTOCOL: the protocol to use for the session url \
  ws -> websocket \
  wss -> websocket secure
- EXTRA_PVPYTHON_ARGS: extra arguments to pass to pvpython (comma-separated, no extra spaces)
  Example: "-dr,--mesa-swr"

Note: If SERVER_NAME and PROTOCOL are not specified, the container defaults to `ws://localhost`.

Paraview & Flask environment variables:

- S3_BUCKET_NAME: The s3 bucket to use for on-the-fly run downloads and for flask server to access for api calls.
- AWS_DEFAULT_REGION: Required by flask for dynamodb access.
- TABLE_NAME: The name of the table that will store run metadata

Note: Each flask route needs only the variables it uses: getTimeSeries needs S3_BUCKET_NAME, availableRuns needs TABLE_NAME, and syncMetadata needs both. A route whose variable is missing returns HTTP 500; /h3lioviz/metadata/health needs none.
The EC2s deployed by swxtrec-cdk set only SERVER_NAME, PROTOCOL and S3_BUCKET_NAME, so availableRuns and syncMetadata fail on them; deployed sites use the swxtrec-cdk Lambda functions for those.
If S3_BUCKET_NAME is not specified, paraview will not be able to download new runs on-the-fly, but will still be able to utilize runs on disk.

### Logging Locations ***Within*** Container

- Flask: `/data/launcher/log/flask.log`
- Paraview: `/data/launcher/log/<session_id>.txt` & `/data/launcher/log/launcherLog.log`
- Apache: `/data/launcher/log/001-pvw_access.log` & `/data/launcher/log/001-pvw_error.log`

## Prerequisites

To work with this repository, you will need:

- [Docker](https://docs.docker.com/get-docker/).
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) (for AWS deployments)
- [Node.js](https://nodejs.org/) 20 (only for building the frontend yourself, steps 1-7 below; the frontend's `package.json` requires `>=20.19.3 <21`)
- Access to the SWx-TREC AWS ECR repository (for deployments). swxtrec-cdk always pulls the image from this repository, in every account; only the image tag is configurable.

## Building the Docker Image

> [!NOTE]
> The Dockerfile now includes steps to build the frontend. You can build the image with an included environment file by selecting one of the currently supported options (`dev`, `prod`, or `swpc`) via a Docker build argument (default is `dev` if unspecified). Other environments (for example, a future `noaa` option) are not yet available. The `--build-arg BUILD_FRONTEND` is optional since this argument defaults to true. 

On an Apple Silicon Mac, add `--platform linux/amd64` to every `docker build` and `docker run` command below: the image downloads x86_64 builds of ParaView, the AWS CLI and s5cmd.

Build the image locally:

```bash
docker build --build-arg FRONTEND_ENVIRONMENT=prod --build-arg BUILD_FRONTEND=true -t h3lioviz .
```

For advanced usage and to build a different version of the frontend locally, you can
follow steps 1-7 to build and modify the frontend code directly.

### Step 1: Clone/Download the Frontend

The frontend is in our [WEBAPPS Bitbucket](https://bitbucket.lasp.colorado.edu/projects/WEBAPPS) for LASP internal users, but there is also a [public mirror available on Github](https://github.com/SWxTREC/h3lioviz).

### Step 2: Update Environment Variables to Point to Your Backend

Enter the repo and open `src/environments/environment.dev.ts` and `src/environments/environment.prod.ts`.

Edit the following configuration items:

- Set `environment.aws.api` to `https://h3lioviz-api.{your-domain}/` (a CloudFront distribution created by swxtrec-cdk, which forwards to a load balancer that routes to AWS Lambda functions).
- Set `environmentConfig.sessionManagerURL` to `https://paraview-web.{your-domain}/h3lioviz/paraview/` (the trailing slash is required; the launcher returns 404 without it).

### Step 3: Install the Frontend Dependencies

```bash
npm install
```

### Step 4: Build the Desired Frontend

Depending on whether or not you're trying to build the frontend for dev or prod run:

```bash
npm run build:dev
```

or

```bash
npm run build:prod
```

### Step 5: Clone the h3lioviz-server Repository

```bash
git clone https://github.com/SWxTREC/h3lioviz-server.git
cd h3lioviz-server
```

### Step 6: Copy the Frontend into h3lioviz-server

Copy everything in the frontends dist directory into pvw/www

Note: After the copy you should have the following directory pvw/www/h3lioviz/

```bash
cp -r {frontend-path}/dist/* pvw/www
```

### Step 7: Build the Docker Image

Build the image locally using a pre-built frontend:

```bash
docker build --build-arg BUILD_FRONTEND=false -t h3lioviz .
```

## Deployment

### Local Development Deployment

This will run the h3lioviz docker image you have build locally. If you want to pull the latest dev/prod image, replace h3lioviz:latest with public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:<tag> with the dev or prod tag.

NOTE: This does not currently entirely work due to the container's dependency on AWS resources.
The paraview code will serve the websocket just fine, and the frontend will be accessible, but the flask routes served by the container as well as on-the-fly run downloading will not work.

```bash
docker run -p 0.0.0.0:8080:80 \
  -e SERVER_NAME=127.0.0.1:8080/h3lioviz \
  -e PROTOCOL=ws \
  -v ${PWD}/test-data:/data \
  -it h3lioviz:latest
```

The server will be available at `http://127.0.0.1:8080`.

To test changes to the pvw directory without rebuilding, add `-v ${PWD}/pvw:/pvw` to mount it directly within the docker container.
The mount hides the image's `/pvw/venv`, so flask will not start, and the frontend bundled in `/pvw/www`, so copy a built frontend into `pvw/www` first (see step 6).
Changes to `pvw/server` take effect in the next ParaView session; changes to `pvw/flask` or `pvw/launcher/config.json` need a container restart.

### AWS ECR Deployment

The official ECR repository is: `public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa`

#### Pushing to Development

1. Authenticate with AWS ECR:

```bash
aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/swx-trec/
```

2. Tag your image:

```bash
docker tag h3lioviz:latest public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev
```

3. Push to ECR:

```bash
docker push public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev
```

Note that the EC2 instance will not update the docker image until it has fully rebooted.
If you want to manually force the update, remotely connect to the instance and run:

```bash
sudo su
export PATH="$PATH:/usr/local/bin"
/docker/docker-launch.sh
```

#### Pushing to Production

Follow the same steps as above, but use the `:prod` tag instead of `:dev`:

```bash
docker tag h3lioviz:latest public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:prod
docker push public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:prod
```

#### Pushing to Legacy (Testing Environment)

For testing on legacy, create an ECR (or use a pre-existing one) and run the above commands referencing that ECR.

1. Authenticate with the legacy account:

```bash
aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin <ECR Name>
```

2. Tag and push:

```bash
docker tag h3lioviz:latest <ECR Name>
docker push <ECR Name>
docker logout public.ecr.aws
```

> **Note**: Remember to log out of the public ECR after pushing, as credentials can interfere with accessing public AWS ECR repositories.

## Test Data

### Downloading Test Data

You can download test data from one of the dev account data buckets. Thinned runs are around 1.2GB.

1. Locate a bucket called `h3lioviz.<domain_name>` (for example `h3lioviz.dev.swx-trec.com`)
2. Identify a run to download: `/data/h3lioviz/pv-ready-data-<run_id>`
3. Download the run: `s5cmd cp "s3://<bucket_name>/data/h3lioviz/pv-ready-data-<run_id>/*" ./test-data/pv-ready-data-<run_id>/`

### Generating Your Own Data

Reference the [README](scripts/README.md) in `scripts/` for instructions on how to generate these runs and make them visible to h3lioviz-server.

### Data Format

The primary data files used by ParaView for visualization are:

- `pv-tim.XXXX.nc` - Time-step files (where XXXX is the time-step number)
- `data_*.vts` - Time-step files for Euhforia runs
- `pv-data-3d.nc` - A single file holding every time step, for runs processed the legacy way

These files contain the 3D heliospheric simulation data.
A run directory also holds the satellite evolution files (`evo.*.json`, which ParaView reads for satellite positions and the frontend graphs; `process_output.py` also writes `evo.*.nc` copies) and `metadata.json`, which describes the run.

## Adding Python Packages

Python packages are installed via a virtual environment during the Docker build process.

1. Add the required packages to `pvw/requirements.txt`
2. Rebuild the Docker image

The virtual environment is created by the Dockerfile and is utilized for flask in server.sh. Note that there is code in place to allow paraview to access the venv, but is currently not in use.

> **Note**: While paraview can get access to additional python packages, our current version of paraview was built with a now outdated SSL version, which makes it close to impossible to use libraries like boto3.
> We currently rely on the s5cmd command line tool for downloading runs. This restriction only applies to python scripts invoked using the pvpython command.

## Repository Structure

```
h3lioviz-server/
├── docker/
│   ├── config/
│   │   └── apache/        # Apache configuration
│   └── scripts/           # Container initialization script
├── docs/                  # More documentation
├── pvw/
│   ├── flask/             # Flask metadata API server
│   ├── launcher/          # ParaView Web launcher configuration
│   ├── server/            # ParaView Web visualization server
│   ├── www/               # Frontend web application
│   └── requirements.txt   # Python dependencies
├── scripts/               # Data processing scripts
├── test-data/             # Test simulation data
├── Dockerfile             # Container build definition
├── LICENSE                # License
└── README.md              # Root documentation page
```

## Container Implementation Details

The container's entrypoint is `/opt/paraviewweb/scripts/server.sh`, which:

1. Updates the paraview-web launcher config (pvw/launcher/config.json) based on docker environment variables
2. Starts flask webserver
3. Starts/restarts apache service
4. Starts the paraview websockets launcher
