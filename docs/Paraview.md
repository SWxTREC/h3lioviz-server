# Paraview

## Paraview Web

### Paraview Web Launcher

The paraview web launcher is what is responsible for starting a new process for each user, generating a session key,
and distributing that session key back to the user and to the visualization server. It listens on localhost:9000. 

#### Configuration

The configuration file is located in [../pvw/launcher/config.json](../pvw/launcher/config.json). Here are a few of the most important parameters:
- configuration.sessionURL: This is the url that the paraview launcher expects requests for new sessions to be made to. SESSION_URL_ROOT is replaced with $PROTOCOL://$SERVER_NAME which are both environment variables set by the docker-compose file mentioned in [../README.md](../README.md).
- configuration.proxyFile: This is the file that will store mappings between sessionIDs and ports that have been opened for the visualization. 
- apps.visualizer: This will run pvpython and is forked for each invocation of the paraview web launcher
- apps.visualizer-mpi: This provides parallelization for batch operations. I am not sure if this is being used. 

#### Logging

`/data/launcher/log/launcherLog.log`

### Paraview Web Visualizer

This is what handles creating the visualizations displayed on the frontend. 
These visualizations are controlled via RPC calls in [../pvw/server/app.py](../pvw/server/app.py). 
Some of the features that come from paraview are enabled and configured by [../pvw/server/app_server.py](../pvw/server/app_server.py) 
along with other parameters necessary to get the server running. 
The features I mentioned are enabled using the [web.protocols module of pvpython](https://www.paraview.org/paraview-docs/latest/python/paraview.web.html).
An example of one such protocol is the ParaviewWebMouseHandler which automatically handles RPC calls for manipulating the visualization on the frontend.

#### Logging

`/data/launcher/log/<hashed_session_id>.log`