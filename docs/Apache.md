# Apache

## Configuration File

The configuration file for the apache server is located in [../docker/config/apache/001-pvw.conf](../docker/config/apache/001-pvw.conf).

## Configuration Details

The Apache server listens on port 80 and proxies the Paraview Launcher, Paraview Web Server, and Flask Server. Here are the rules for those proxies:

- Paraview Launcher Proxy: localhost:80/h3lioviz/paraview -> http://localhost:9000/paraview
- Paraview Web Server Proxy: localhost:80/h3lioviz/proxy?sessionId=XXXX&path=ws -> ws://SOME_PORT:SESSION_ID/PATH
- Flask Server Proxy: localhost:80/h3lioviz/metadata -> http://localhost:5000/h3lioviz/metadata

