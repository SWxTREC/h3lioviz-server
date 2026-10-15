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
- `/h3lioviz/paraview/` -> ParaView Web service (creates visualization sessions; the trailing slash is required)
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

Note: If either SERVER_NAME or PROTOCOL is unset, the container uses `ws://localhost`.

Paraview & Flask environment variables:

- S3_BUCKET_NAME: The s3 bucket to use for on-the-fly run downloads and for flask server to access for api calls.
- AWS_DEFAULT_REGION: Required by flask for dynamodb access.
- TABLE_NAME: The name of the table that will store run metadata

Note: Each flask route needs only the variables it uses: getTimeSeries needs S3_BUCKET_NAME, availableRuns needs TABLE_NAME and AWS_DEFAULT_REGION, and syncMetadata needs all three. A route whose variable is missing returns HTTP 500; /h3lioviz/metadata/health needs none.
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
- [Node.js](https://nodejs.org/) 20 (only for building the frontend yourself, steps 1-7 below; the frontend's `package.json` declares `>=20.19.3 <21`; npm only warns outside that range, and the Dockerfile builds the frontend with Node 22)
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

This will run the h3lioviz docker image you have built locally. To run a published image instead, replace h3lioviz:latest with public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:<tag>. The `dev` and `prod` images carry no frontend of their own (deployed sites get it from CloudFront). With the `pvw` mount, the server code and any frontend come from your checkout, whichever image you run.

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

The official ECR repository is: `public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa`. It is a public repository in the dev account. Each EC2 runs `/docker/docker-launch.sh` at boot, which pulls the tag its account names in swxtrec-cdk's `cdk.json` (`paraview_image_tag`: `prod` for prod and swpc; `dev` for dev, legacy and any account that leaves it unset). A new image therefore goes live at the instance's next start.

#### Development and production: GitHub Actions

`dev` and `prod` are built and moved by workflows, not by hand:

- **Build and publish image** (`.github/workflows/build-image.yml`) runs on every push to `main` that changes the image (the `Dockerfile`, `pvw/` or `docker/`) or the build workflow and smoke test; a docs-only merge builds nothing. It builds the image without a bundled frontend (`BUILD_FRONTEND=false`; deployed sites get the frontend from CloudFront), checks that it starts a ParaView session and that its metadata server answers the health route (`.github/scripts/smoke-test.sh`), pushes it under a tag no other build reuses (`sha-<12 hex of the commit>`, or `sha-<12 hex>-run<run id>.<attempt>` for a rebuild of the same commit), and moves `dev` to it unless `main` has taken a newer change to the image in the meantime. Pull requests run the build and the smoke test only. A manual dispatch on `main` from the Actions tab always builds, and publishes as a push does.
- **Promote image to prod** (`.github/workflows/promote-image.yml`) is started by hand from the Actions tab once an image has been tested on dev, and waits for a reviewer's approval (the `prod` environment; see `.github/iam/README.md`). Give it that image's digest (`sha256:...`, shown in the build's summary) or build tag; an empty input promotes whatever `dev` names at that moment, which is a newer image if main has moved since the test. Nothing is rebuilt, so prod runs exactly the image promoted. Rolling back is the same dispatch with an older digest.

Older images are kept for rollback and reference: each build keeps its own tag, and each promotion that moves `prod` tags the image it named before as `prod-until-<UTC time>-run<run id>`. The workflows never delete images. Very old ones can be removed by hand; see `.github/iam/README.md` for which must stay. The workflows' AWS roles are also described there.

To run a new image before the instance's next start, remotely connect to the instance and run:

```bash
sudo su
export PATH="$PATH:/usr/local/bin"
/docker/docker-launch.sh
```

#### Pushing by hand

Other tags (for example the SWPC stand-in images `noaa` and `noaatest`, which bundle their own frontend build) are still pushed by hand. Do not push `dev`, `prod`, `sha-*` or `prod-until-*` by hand: the workflows own those tags, a hand push records no commit, overwriting a `sha-` tag would make a promotion by that tag pick up a different image, and overwriting a `prod-until-` tag would lose track of an earlier prod image.

1. Authenticate with AWS ECR:

```bash
aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws
```

2. Tag your image:

```bash
docker tag h3lioviz:latest public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:<tag>
```

3. Push to ECR:

```bash
docker push public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:<tag>
```

#### Testing an image on legacy or a test account

Legacy and any account that leaves `paraview_image_tag` unset pull `dev` from the swx-trec repository, so they follow the build workflow. To run a different image there, push it to swx-trec under a tag of your own (as above), then connect to the instance (Session Manager is enabled on it), set the `image:` line in `/docker/docker-compose-cpu.yaml` to your tag, and run `/docker/docker-launch.sh` as shown above, which pulls the image and restarts the container. Set the line back the same way when you are done. Changing `paraview_image_tag` in swxtrec-cdk's `cdk.json` alone does not do this on an existing instance: the tag is written into the instance's user data, which runs only at first boot, so a deploy restarts the instance without rewriting the compose file. It takes effect on a newly created instance.

> **Note**: Log out of the public ECR after pushing (`docker logout public.ecr.aws`), as credentials can interfere with pulling from public AWS ECR repositories.

## Test Data

### Downloading Test Data

You can download test data from one of the dev account data buckets. Thinned runs are around 1.2GB.

1. Locate a bucket called `h3lioviz.<domain_name>` (for example `h3lioviz.dev.swx-trec.com`)
2. Identify a run to download: `data/h3lioviz/pv-ready-data-<run_id>`
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

The virtual environment `/pvw/venv` is created by the Dockerfile and is used for flask in server.sh. ParaView's `pvpython` runs with a separate venv, `/opt/pv-venv` (`pvpython --venv`), which holds wslink and its dependencies, because ParaView 6 no longer bundles it.

> **Note**: Python scripts run with `pvpython` download runs with s5cmd (`pvw/server/app.py`) rather than boto3. ParaView 5.10, which this server ran before the 6.1 upgrade, was built with an outdated SSL version that made libraries like boto3 close to impossible to use there; this has not been re-checked on 6.1.

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
