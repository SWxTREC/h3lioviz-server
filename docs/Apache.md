# Apache

## Configuration File

The configuration file for the apache server is located in [../docker/config/apache/001-pvw.conf](../docker/config/apache/001-pvw.conf).

## Configuration Details

The Apache server listens on port 80 and acts as a reverse proxy for the Paraview Launcher, Paraview Web Servers, and Flask Server.
The frontend root directory is at `../pvw/www`. Apache is also configured to attach the Access-Control-Allow-Origin "\*" header to all traffic it routes to avoid any CORS issues.
Apache determines which port to route paraview web server requests to using /opt/launcher/proxy-mapping.txt

### Proxy Mappings and Routes:

#### Proxy Mappings:

- Paraview Launcher Reverse Proxy: http://localhost:80/h3lioviz/paraview -> http://localhost:9000/paraview
- Paraview Web Servers Reverse Proxy: wss://localhost:80/h3lioviz/proxy?sessionId=XXXX&path=ws -> ws://0.0.0.0:<SOME_PORT>/ws \
  Note that <SOME_PORT> is determined by the sessionId.
- Flask Server Reverse Proxy: http://localhost:80/h3lioviz/metadata -> http://localhost:5000/h3lioviz/metadata

#### Routes and Redirects

- / ->(302) /h3lioviz
- /helioviz ->(302) /h3lioviz
- /h3lioviz/[FRONTEND-FILENAME]
- /h3lioviz/[FRONTEND-DIR]/[FRONTEND-FILENAME]
- If none of the above conditions are included and the request does not hit a reverse proxy endpoint, serve /h3lioviz/index.html.

## Logging

The access log contains the following information:

- Basic information on all incoming http requests to port 80 on the instance including IP address, request method, date, browser, and status code
- It will NOT contain information on web sockets connections unless `LogLevel error proxy:trace5` is added to the config. This can create huge log files, so it is currently off by default.

`/var/log/apache2/001-pvw_access.log`

The error log contains the following information:

- Failed proxies where the downstream server can't be accessed
- Other Apache erorrs

`/var/log/apache2/001-pvw_error.log`
