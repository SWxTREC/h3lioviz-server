# Paraview

## Paraview Web

### Paraview Web Launcher

The paraview web launcher is responsible for starting a new process for each user, generating a session key,
distributing that session key back to the client, and mapping session key to the port the forked visualization server is running on. It listens on localhost:9000.

#### Configuration

The configuration file is located in [../pvw/launcher/config.json](../pvw/launcher/config.json). Here are a few of the most important parameters:

- configuration.sessionURL: This is the url that the paraview launcher expects requests for new sessions to be made to. SESSION_URL_ROOT is replaced with $PROTOCOL://$SERVER_NAME by server.sh. Both $PROTOCOL and $SERVER_NAME are environment variables set by the docker-compose file mentioned in [../README.md](../README.md).
- configuration.proxyFile: This is the file that will store mappings between sessionIDs and ports that have been opened for the websocket controlled visualization.
- apps.visualizer: This will run pvpython and is forked for each invocation of the paraview web launcher
- apps.visualizer-mpi: This provides parallelization for batch operations. I am not sure if this is being used.

#### Logging

`/data/launcher/log/launcherLog.log`

### Paraview Web Visualizer/Server

Each connection will create a new paraview web server process with it's own port.
The capablities of the visualizer are determined by which VtkWebProtocols are registered to the paraview server.
An example of one such protocol is the ParaviewWebMouseHandler which automatically handles RPC calls for moving the visualization on the frontend with your mouse.
Another example is the custom ParaviewWebProtocol, a subclass of vtkWebProtocol, defined in [../pvw/server/app.py](../pvw/server/app.py).
This is how we define custom actions on the visualization that the frontend can request via RPC calls.
The VtkWebProtocols are registered with the paraview web server in [../pvw/server/app_server.py](../pvw/server/app_server.py).
The app_server.py file also configures some of these protocols where necessary.

#### Logging

`/data/launcher/log/<hashed_session_id>.log`

## Paraview GUI

### Installation

Navigate to the [Paraview download page](https://www.paraview.org/download/). Set the version selector to v5.10 and change your platform appropriately (the platforms are below the version selector).

### Basic Usage

You will need processed data in order to use the Paraview GUI for visualizing h3lioviz data. There are instructions for this in [../scripts/README.md](../scripts/README.md). It is recomended you use the default downsampling as the visualization can be pretty hard on your system depending on it's specs.

Once you have processed data open the Paraview GUI and select the 'Open File' icon. Then navigate to your pv-ready-{RUN_ID} folder and select the pv-tim-\*.nc files. You should be albe to select them all at once. Once you open these files you will be prompted for a reader. Select 'NetCDF Reader'. This will read the data into Paraview. In the left pannel navigate to 'Properties' and select apply.

Right click on the pv-tim.0\* files in the 'Pipeline Browser' and select 'Add Filter' then 'Alphabetical' and select 'Cell Data to Point Data'. Next navigate to the 'Properties' tab and select 'Apply'. This filter will convert from blocky cells of data to a smoother representation.

After applying the Cell Data to Point Data filter go to the bar above that panel and select 'Slice'. Navigate to 'Properties' if it isn't open already. Hit the 'Z Normal' button under 'Plane Parameters' in order to snap the slice to the correct plane. Then hit 'Apply'.

In order to see the data you will need to select a field. You can do this by editing the 'Display' settings. Under 'Coloring' you can change the 'Solid color' to be any variable you want from the data. Then hit 'Apply'

Finally you can play back the files using the controls in the top bar.
