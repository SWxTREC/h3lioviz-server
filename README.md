# Enlil 3D Server

This repository hosts the code to run a Paraview server for
the 3D Enlil output.

## Running the server

To run the server locally on port 1234 the command would be something like the following.

```bash
/path/to/paraview/bin/pvpython pvw/server/pv_server.py --port 1234
```

The frontend can then use wslink to use a websocket to connect to port 1234.

## Installation

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
pvpython pvw/server/pv_server_statefile.py --port 1234 --file /path/to/test.nc
```

where the port is `1234` for local development, and the path to the input file is
given to the `--file` command line argument.

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
    docker run -p 0.0.0.0:9000:80 -e SERVER_NAME="127.0.0.1:9000" -e PROTOCOL="ws" -v ${PWD}/pvw:/pvw -v ${PWD}/data:/data /path/to/frontend/dist/swt:/frontend -it pvw-enlil-osmesa
    ```

    The container requires several input volumes that contain the frontend
    code mounted at `/frontend` inside the container, a data directory that contains
    a `test.nc` input file mounted to the `/data` location inside the
    container, and the `pvw` directory from this repository mounted
    at the `/pvw` location.
