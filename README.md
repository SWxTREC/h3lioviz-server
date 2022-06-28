# H3LIO Viz

This repository hosts the code to run a Paraview server for
the 3D heliospheric output, from codes such as Enlil and EUHFORIA.

## Local Development and Testing

To test the code locally, you will need the Docker container that has Paraviewweb installed, the test data files, and this cloned git repository. This will require
Docker and the aws-cli programs to run locally.

[Get Docker](https://docs.docker.com/get-docker/)
[Get aws-cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)

### Getting the Docker container locally

Use your AWS credentials to authorize yourself to the AWS Registry.

```bash
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 135080795405.dkr.ecr.us-east-1.amazonaws.com
```

Pull the container to your local system. (Note: the `osmesa` container is for CPU rendering, if you want to test GPU rendering you'll need the `egl` container)

```bash
docker pull 135080795405.dkr.ecr.us-east-1.amazonaws.com/enlil:pvw-enlil-osmesa
```

### Getting the test data locally

This is a 4 GB archive containing the relevant data products. You will need AWS credentials to access this bucket as well, which can be given as a profile `--profile swx-trec-base`

```bash
aws s3 cp s3://enlil-data/pv-ready-data.tar .
tar -xf pv-ready-data.tar
```

### Running the package locally

```bash
docker run -p 0.0.0.0:8080:80 -e SERVER_NAME=127.0.0.1:8080 -e PROTOCOL=ws -v ${PWD}/pvw:/pvw -v ${PWD}/pv-ready-data:/data -it 135080795405.dkr.ecr.us-east-1.amazonaws.com/enlil:pvw-enlil-osmesa
```

## Building the Dockerfile

This was built following the [guides](https://github.com/Kitware/paraviewweb/tree/master/tools/docker) from Kitware and their repositories. The `dockerfile/` directory contains all of the content that needs to be
copied into the container.

> **_NOTE:_**  You may need to make the scripts executable if they weren't cloned with those permissions.
> `chmod a+x scripts/`

1. Download the Paraview binary for linux and save it into the `docker/binaries/` directory. [https://www.paraview.org/download/](https://www.paraview.org/download/)
    > **_NOTE:_**  If you have a GPU available use the `egl` version, otherwise use the `osemsa`.

2. Build the image locally with an appropriate tag.

    ```bash
    cd docker
    docker build --rm -t pvw-enlil-osmesa .
    ```

3. Run the image setting the proper environment variables and mounting the proper directories.

    ```bash
    docker-compose --env-file docker/.env.local --file docker/docker-compose.yaml up
    ```

    The container requires several input volumes that contain the data directory (that contains
    a `pv-data-3d.nc` input file) mounted at the `/data` location inside the
    container, and the `pvw` directory from this repository mounted
    at the `/pvw` location. To handle that, you can create a `.env.local` environment variable
    file to point to these data locations on your local system. For instance, your local file may
    look like:

    ```bash
    PVW_BACKEND=/home/code/enlil-3d-server/pvw  
    PVW_DATA=/home/data/pv-ready-data
    ```

## Running the server without Docker

To run the server locally on port 1234 the command would be something like the following.

```bash
/path/to/paraview/bin/pvpython pvw/server/pv_server.py --port 1234
```

The frontend can then use wslink to use a websocket to connect to port 1234.

### Installation

The standard Paraview binaries seem to have issues with missing packages. Creating a
new virtual environment and installing the dependencies manually seemed to be the easiest
way to address that. The below commands worked for me to get a Paraview 5.9.1 server running.
The conda pvpython binary is missing `autobahn` and `wslink`, the latter of which isn't
in the conda package manager, so I had to use pip to install it.

```bash
conda create -n paraview python=3.9 paraview autobahn
pip install wslink
```

Now, you can run the server with

```bash
pvpython pvw/server/app_server.py --port 1234 --dir /path/to/<pv-data-3d.nc>
```

where the port is `1234` for local development, and the path to the input file is
given to the `--dir` command line argument.
