# H3lioViz

This repository hosts the code to run a ParaView visualization server for 3D heliospheric output from codes such as Enlil and Euhforia, as well as a way to fetch the associated metadata and host the frontend.

## Overview

The H3lioViz server is a containerized application that provides:
- **ParaView Web Service**: 3D visualization of heliospheric simulation data
- **Flask Metadata API**: Endpoints for retrieving run information and time-series data
- **Apache Web Server**: Serves the frontend and acts as a reverse proxy

### Architecture

Apache serves as the main entry point and handles routing:
- `/` and `/h3lioviz/` → Frontend web application
- `/paraview` → ParaView Web service (creates visualization sessions)
- `/proxy` → WebSocket proxy (maps session IDs to ParaView ports)
- `/metadata/` → Flask API for metadata operations

### Flask API Endpoints

The Flask server provides the following REST endpoints:

- `GET /metadata/health` - Health check endpoint, returns `{"status": "ok"}`
- `GET /metadata/getTimeSeries/<run_id>/<satellite>` - Retrieves time-series data for a specific run and satellite
- `GET /metadata/availableRuns` - Lists all available simulation runs
- `GET /metadata/syncMetadata` - Synchronizes metadata from data sources

## Prerequisites

To work with this repository, you will need:

- [Docker](https://docs.docker.com/get-docker/).
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) (for AWS deployments)
- Access to the SWx-TREC AWS ECR repository (for production/development deployments), or a personal one (for deployments to an account other than prod or dev).

## Building the Docker Image

### Step 1: Clone the Repository

```bash
git clone https://github.com/SWxTREC/h3lioviz-server.git
cd h3lioviz-server
```

### Step 2: Download & Extract the ParaView Binaries

Download the ParaView 5.10.1 binaries for Linux and extract the tarball into the `docker/binaries/` directory:

```bash
curl https://www.paraview.org/paraview-downloads/download.php\?submit\=Download\&version\=v5.10\&type\=binary\&os\=Linux\&downloadFile\=ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz --output ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz
tar -xzvf ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz -C ./docker/binaries/ --strip-components=1
```

You should see `bin`, `lib`, and `share` directories in `docker/binaries/`.

> **Note**: Use the `osmesa` version for CPU rendering. As of October 2025, we have not successfully run a GPU-based (egl) version of paraview. 

### Step 3: Build the Docker Image

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
  -e SERVER_NAME=127.0.0.1:8080 \
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
aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev
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

TODO: Document how to download test data

### Data Format

The primary data files used by ParaView for visualization are:
- `pv-tim.XXXX.nc` - Time-step files (where XXXX is the time-step number)

These NetCDF files contain the 3D heliospheric simulation data. The test-data directory contains additional files (satellite evolution files, metadata, etc.) that are graphed by the frontend but are not directly loaded by ParaView.

## Adding Python Packages

Python packages are installed via a virtual environment during the Docker build process.

1. Add the required packages to `pvw/requirements.txt`
2. Rebuild the Docker image

The virtual environment is created and activated automatically by the container's entrypoint script (`docker/scripts/server.sh`).

> **Note**: There is currently code in place to configure the venv to allow for paraview to access it, but it has an outdated embedded SSL version, which makes it close to impossible to use libraries like boto3. We currently rely on using AWS CLI commands for downloading runs.

## Repository Structure

```
h3lioviz-server/
├── docker/
│   ├── binaries/          # ParaView binaries (bin/, lib/, share/)
│   ├── config/
│   │   └── apache/        # Apache configuration
│   └── scripts/           # Container initialization scripts
├── pvw/
│   ├── flask/             # Flask metadata API server
│   ├── launcher/          # ParaView Web launcher configuration
│   ├── server/            # ParaView Web visualization server
│   ├── www/               # Frontend web application
│   └── requirements.txt   # Python dependencies
├── scripts/               # Data processing scripts
├── test-data/             # Test simulation data
├── Dockerfile             # Container build definition
└── docker-compose-local.yaml  # Local testing configuration
```

## Container Implementation Details

> **TODO**: This section will be updated once the container initialization flow and ParaView Python code are finalized. Current implementation details may change.

The container's entrypoint is `/opt/paraviewweb/scripts/server.sh`, which:
1. Sets up the Python virtual environment
2. Configures Apache routing via `addEndpoints.sh`
3. Launches the ParaView Web launcher via `start.sh`

TODO: Document the complete initialization flow, including:
- How the launcher manages ParaView sessions
- Port allocation and session management
- WebSocket proxy configuration
- Container networking and EC2 integration
