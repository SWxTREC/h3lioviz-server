# Enlil 3D Server

This repository hosts the code to run a Paraview server for
the 3D Enlil output.

## Running the server

To run the server locally on port 1234 the command would be something like the following.

```bash
/path/to/paraview/bin/pvpython pv_server_statefile.py --port 1234
```

The frontend can then use wslink to use a websocket to connect to port 1234.
