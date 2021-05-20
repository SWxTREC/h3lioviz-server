# Enlil 3D Server

This repository hosts the code to run a Paraview server for
the 3D Enlil output.

## Running the server

To run the server locally on port 1234 the command would be something like the following.

```bash
/path/to/paraview/bin/pvpython pv_server.py --port 1234
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
pvpython pv_server_statefile.py --port 1234
```
