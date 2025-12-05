# H3lioViz

This repository hosts the code to run a ParaView visualization server for 3D heliospheric output from codes such as Enlil and Euhforia, as well as a way to fetch the associated metadata and host the frontend.

Live deployment for interactive use: <https://swx-trec.com/h3lioviz/>

## Overview

The H3lioViz server is a containerized application that provides:
- **ParaView Web Service**: 3D visualization of heliospheric simulation data
- **Flask Metadata API**: Endpoints for retrieving run information and time-series data
- **Apache Web Server**: Serves the frontend and acts as a reverse proxy

### Architecture

Apache serves as the main entry point and handles routing:
- `/` -> Reroutes to `/h3lioviz/`
- `/h3lioviz/` -> Frontend web application
- `/h3lioviz/paraview` -> ParaView Web service (creates visualization sessions)
- `/h3lioviz/proxy` -> WebSocket proxy (maps session IDs to ParaView ports)
- `/h3lioviz/metadata/` -> Flask API for metadata operations

### Flask API Endpoints

The Flask server provides the following REST endpoints:

- `GET /h3lioviz/metadata/health` - Health check endpoint, returns `{"status": "ok"}`
- `GET /h3lioviz/metadata/getTimeSeries/<run_id>/<satellite>` - Retrieves time-series data for a specific run and satellite
- `GET /h3lioviz/metadata/availableRuns` - Lists all available simulation runs
- `GET /h3lioviz/metadata/syncMetadata` - Synchronizes metadata from data sources

### Container Startup Script
The docker/scripts/server.sh script initializes the docker container by using certain environment variables to update configuration files.

server.sh (paraview/websockets) environment variables:
  - SERVER_NAME: the server name to use for the session url. Gets returned from paraview-web to tell the frontend where to connect to the websocket. Note that the current routing expects {domain}/h3lioviz.
  - PROTOCOL: the protocol to use for the session url 
      ws -> websocket
      wss -> websocket secure
  - EXTRA_PVPYTHON_ARGS: extra arguments to pass to pvpython (comma-separated, no extra spaces)
      Example: "-dr,--mesa-swr"
    
Note: If SERVER_NAME and PROTOCOL are not specified, the container defaults to `ws://localhost`.

Paraview & Flask environment variables:
  - S3_BUCKET_NAME: The s3 bucket to use for on-the-fly run downloads and for flask server to access for api calls.
  - AWS_DEFAULT_REGION: Required by flask for dynamodb access.
  - TABLE_NAME: The name of the table that will store run metadata

Note: If the two above parameters are not specified, any calls to flask (except /h3lioviz/metadata/health) will fail. If S3_BUCKET_NAME is not specified, paraview will not be able to download new runs on-the-fly, but will still be able to utilize runs on disk.

### Logging Locations Within Container
 - Flask: `/data/launcher/log/flask.log`
 - Paraview: `/data/launcher/log/<hashed_session_id>.log` & `/data/launcher/log/launcherLog.log`
 - Apache: `/var/log/apache2/001-pvw_access.log` & `/var/log/apache2/001-pvw_error.log`

## Prerequisites

To work with this repository, you will need:

- [Docker](https://docs.docker.com/get-docker/).
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) (for AWS deployments)
- Access to the SWx-TREC AWS ECR repository (for production/development deployments), or a personal one (for deployments to an account other than prod or dev).

## Building the Docker Image

### Step 1: Clone/Download the Frontend
The frontend is in our [WEBAPPS Bitbucket](https://bitbucket.lasp.colorado.edu/projects/WEBAPPS) for LASP internal users, but there is also a [public mirror available on Github](https://github.com/SWxTREC/h3lioviz).

### Step 2: Update Environment Variables to Point to Your Backend

Enter the repo and open `src/environments/environment.dev.ts` and `src/environemnts/environment.prod.ts`.

Edit environment.aws.api and environment.aws.api to `https://h3lioviz-api.{your-domain}/` and environmentConfig.sessionManagerURL: `https://paraview-web.{your-domain}/paraview`. 

### Step 3: Install the Frontend Dependencies
```bash
npm install
```

### Step 4: Build the Desired Frontend
Depending on whether or not you're trying to build the frontend for dev or prod run: 

```bash
npm build:dev
```

or

```bash
npm build:prod
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

Build the image locally:

```bash
docker build -t h3lioviz .
```

## Deployment

### Local Development Deployment
This will run the h3lioviz docker image you have build locally. If you want to pull the latest dev/prod image, replace h3lioviz:latest with public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:<tag> with the dev or prod tag.

NOTE: This does not currently entirely work due to the container's dependency on AWS resources. The paraview code will serve the websocket just fine, but the frontend served by the container will not function. See [Building & Running h3lioviz-server](https://confluence.lasp.colorado.edu/spaces/MODSDB/pages/271520419/Building+Running+h3lioviz-server) for more details.

```bash
docker run -p 0.0.0.0:8080:80 \
  -e SERVER_NAME=127.0.0.1:8080/h3lioviz \
  -e PROTOCOL=ws \
  -v ${PWD}/pvw:/pvw \
  -v ${PWD}/test-data:/data \
  -it h3lioviz:latest
```

The server will be available at `http://127.0.0.1:8080`.

The pvw directory is mounted directly within the docker container, so any changes to the code within will be reflected next time the connection is refreshed.

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
Note that the EC2 instance will not update the docker image until it has fully rebooted. If you want to manually force the update, remotely connect to the instance and run:
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

> **Note**: Remember to log out of the public ECR after pushing, as credentials can interfere with other AWS Docker processes.

## Test Data

### Downloading Test Data
You can download test data from one of the dev account data buckets. Thinned runs are around 1.2GB.

1. Locate a bucket called `h3lioviz.<domain_name>.com`
2. Identify a run to download: `/data/h3lioviz/pv-ready-data-<run_id>`
3. Download the run: `s5cmd cp "s3://<bucket_name>/data/h3lioviz/pv-ready-data-<run_id>/\*" ./test-data/pv-ready-data-<run_id>/`

### Generating Your Own Data

Reference the [README](scripts/README.md) in `scripts/` for instructions on how to generate these runs and make them visible to h3lioviz-server. 

### Data Format

The primary data files used by ParaView for visualization are:
- `pv-tim.XXXX.nc` - Time-step files (where XXXX is the time-step number)

These NetCDF files contain the 3D heliospheric simulation data. The test-data directory contains additional files (satellite evolution files, metadata, etc.) that are graphed by the frontend but are not directly loaded by ParaView.

## Adding Python Packages

Python packages are installed via a virtual environment during the Docker build process.

1. Add the required packages to `pvw/requirements.txt`
2. Rebuild the Docker image

The virtual environment is created by the Dockerfile and is utilized for flask in server.sh. Note that there is code in place to allow paraview to access the venv, but is currently not in use.

> **Note**: While paraview can get access to additional python packages, our current version of paraview was built with a now outdated SSL version, which makes it close to impossible to use libraries like boto3. We currently rely on using AWS CLI commands for downloading runs. This restriction only applies to python scripts invoked using the pvpython command.

## Repository Structure

```
h3lioviz-server/
├── docker/
│   ├── binaries/          # ParaView binaries (bin/, lib/, share/)
│   ├── config/
│   │   └── apache/        # Apache configuration
│   └── scripts/           # Container initialization script
├── pvw/
│   ├── flask/             # Flask metadata API server
│   ├── launcher/          # ParaView Web launcher configuration
│   ├── server/            # ParaView Web visualization server
│   ├── www/               # Frontend web application
│   └── requirements.txt   # Python dependencies
├── scripts/               # Data processing scripts
├── test-data/             # Test simulation data
├── Dockerfile             # Container build definition
```

## Container Implementation Details

> **TODO**: This section will be updated once the container initialization flow and ParaView Python code are finalized. Current implementation details may change.

The container's entrypoint is `/opt/paraviewweb/scripts/server.sh`, which:
1. Updates the paraview-web launcher config (pvw/launcher/config.json) based on docker environment variables
2. Starts flask webserver
3. Starts/restarts apache service
4. Starts paraview-web service.



TODO: Document paraview-web execution details:
- How the launcher manages ParaView sessions
- Port allocation and session management
