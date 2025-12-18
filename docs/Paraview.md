# Paraview

## Paraview Web

### Paraview Web Launcher

The paraview web launcher is responsible for starting a new process for each user, generating a session key,
distributing that session key back to the client, and mapping session key to the port the forked visualization server is running on.
It listens on localhost:9000. Note that the launcher is only communicated to via HTTP. Our cloud infrastructure handles any HTTP-to-HTTPS
encapsulation. The Paraview web server/visualizer is what handles all websocket connections.

#### Configuration

The configuration file is located in [../pvw/launcher/config.json](../pvw/launcher/config.json). Here are a few of the most important parameters:

- configuration.sessionURL: This is the url that the paraview launcher expects requests for new sessions to be made to. SESSION_URL_ROOT is replaced with \$PROTOCOL://\$SERVER_NAME by server.sh. Both $PROTOCOL and $SERVER_NAME are environment variables that need to be set by a docker-compose file or command line flags.
- configuration.proxyFile: This is the file that will store mappings between sessionIDs and ports that have been opened for the websocket controlled visualization.
- apps.visualizer: This will run pvpython and is forked for each invocation of the paraview web launcher
- apps.visualizer-mpi: This provides parallelization for batch operations. I am not sure if this is being used.

#### Logging

The following logs contain:

- Launcher requests for a new session via POST to /h3lioviz/paraview
- Status requests which the frontend uses to see if the server is online. These are GET requests to /h3lioviz/paraview?{some_random_number}.
  A 400 status_code indicates the server is accessible
- Any errors/warnings during the session launching process

`/data/launcher/log/launcherLog.log`

### Paraview Web Visualizer/Server

Each connection will create a new paraview web server process with it's own port.
The capabilities of the visualizer are determined by which VtkWebProtocols are registered to the paraview server.
An example of one such protocol is the ParaviewWebMouseHandler which automatically handles RPC calls for moving the visualization on the frontend with your mouse.
Another example is the custom ParaviewWebProtocol, a subclass of vtkWebProtocol, defined in [../pvw/server/app.py](../pvw/server/app.py).
This is how we define custom actions on the visualization that the frontend can request via RPC calls.
The VtkWebProtocols are registered with the paraview web server in [../pvw/server/app_server.py](../pvw/server/app_server.py).
The app_server.py file also configures some of these protocols where necessary.

#### Logging

The following logs contain:

- Rendering logs including the time each step of the render took
- Errors that occurred during the rendering or any errors that resulted from websocket RPC calls to the server
- Any other information output by the python files in `../pvw/server` including information on on the fly run downloads

`/data/launcher/log/<hashed_session_id>.log`

#### Testing Websockets Without a Frontend

[Install Postman](https://www.postman.com/downloads/) or another API testing tool. You can avoid all their sign in stuff if you don't want to persist your session. 
The value of the id key doesn't seem to matter much, but it must be filled out in order for requests to be successful (flingus:0 seems to work just fine).

1. Send a POST request to `https://paraview-web.{YOUR_SUB_DOMAIN}.swx-trec.com/h3lioviz/paraview` with the following JSON body:
```json
{
    "application": "visualizer",
    "sessionManagerURL": "https://paraview-web.noaa-demo.swx-trec.com/h3lioviz/paraview/"
}
```
Response:
```json
{
    "sessionManagerURL": "https://paraview-web.noaa-demo.swx-trec.com/h3lioviz/paraview/",
    "id": "d6a03303-d6c6-11f0-b000-e555d37416b1",
    "host": "0.0.0.0",
    "port": 9100,
    "sessionURL": "wss://paraview-web.noaa-demo.swx-trec.com/h3lioviz/proxy?sessionId=d6a03303-d6c6-11f0-b000-e555d37416b1&path=ws"
}
```
2. Create a websockets request. Input the url from sessionURL then hit "Connect". You are now connected to the websocket!
3. Send the following test message:
```json
 {"wslink": "1.0", "id": "system:c0:0", "method": "wslink.hello", "args": [{"secret": "wslink-secret"}]}
 ```
Response:
```json
{
    "wslink": "1.0",
    "id": "system:c0:0",
    "result": {
        "clientID": "c387976380d7746a7a42428f388e05ce1"
    }
}
```
4. Before you can do most other RPC calls you need to load a model. Do so with the following message:
```json
{"wslink": "1.0", "id": "rpc:c0ffd038671824f0db5efe17a02f5d3d4:9", "method": "pv.h3lioviz.load_model", "args": ["57758"]}
```
Response:
```json
{
    "wslink": "1.0",
    "id": "rpc:c0ffd038671824f0db5efe17a02f5d3d4:9",
    "result": null
}
```
5. Make any other RPC calls you want. Note that the order of arguments is the same as the function signature below the rpc decorator:
{"wslink": "1.0", "method": "pv.h3lioviz.visibility", "id": "rpc:c0ffd038671824f0db5efe17a02f5d3d4:1", "args": ["color_bar", "off"]}


## Paraview GUI

### Installation

Navigate to the [Paraview download page](https://www.paraview.org/download/). Set the version selector to v5.10 and change your platform appropriately (the platforms are below the version selector).

### Basic Usage

You will need processed data in order to use the Paraview GUI for visualizing h3lioviz data. There are instructions for this in [../scripts/README.md](../scripts/README.md).
It is recommended you use the default downsampling as the visualization can be pretty hard on your system depending on it's specs.

Once you have processed data open the Paraview GUI and select the 'Open File' icon. Then navigate to your pv-ready-{RUN_ID} folder and select the pv-tim-\*.nc files.
You should be able to select them all at once.
Once you open these files you will be prompted for a reader. Select 'NetCDF Reader'. This will read the data into Paraview.
In the left panel navigate to 'Properties' and select apply.

Right click on the pv-tim.0\* files in the 'Pipeline Browser' and select 'Add Filter' then 'Alphabetical' and select 'Cell Data to Point Data'.
Next navigate to the 'Properties' tab and select 'Apply'. This filter will convert from blocky cells of data to a smoother representation.

After applying the Cell Data to Point Data filter go to the bar above that panel and select 'Slice'.
Navigate to 'Properties' if it isn't open already. Hit the 'Z Normal' button under 'Plane Parameters' in order to snap the slice to the correct plane. Then hit 'Apply'.

In order to see the data you will need to select a field. You can do this by editing the 'Display' settings. Under 'Coloring' you can change the 'Solid color' to be any variable you want from the data. Then hit 'Apply'

Finally you can play back the files using the controls in the top bar.
