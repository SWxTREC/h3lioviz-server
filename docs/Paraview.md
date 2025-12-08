# Paraview

## Paraview Web

### Paraview Web Launcher

The paraview web launcher is what is responsible for starting a new processs for each user, generating a session key,
and distributing that session key back to the user and to the visualization server.

### Paraview Web Visualizer

This is what handles creating the visualizations displayed on the frontend. 
These visualizations are controlled via RPC calls in [../pvw/server/app.py](../pvw/server/app.py). 
Some of the features that come from paraview are enabled and configured by [../pvw/server/app_server.py](../pvw/server/app_server.py) 
along with other parameters necessary to get the server running. 
The features I mentioned are enabled using the [web.protocols module of pvpython](https://www.paraview.org/paraview-docs/latest/python/paraview.web.html).
An example of one such protocol is the ParaviewWebMouseHandler which automatically handles RPC calls for manipulating the visualization on the frontend.

## Paraview GUI

