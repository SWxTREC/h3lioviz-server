# Apache

## Configuration File

The configuration file for the apache server is located in [../docker/config/apache/001-pvw.conf](../docker/config/apache/001-pvw.conf).

## Configuration Details

The Apache server listens on port 80 and proxies the Paraview Launcher, Paraview Web Server, and Flask Server. It is also responsible for serving the frontend. Apache attaches Access-Control-Allow-Origin "*" to all traffic it routes. 

### Proxy Mappings:

- Paraview Launcher Proxy: localhost:80/h3lioviz/paraview -> http://localhost:9000/paraview
- Paraview Web Server Proxy: localhost:80/h3lioviz/proxy?sessionId=XXXX&path=ws -> ws://SOME_PORT:SESSION_ID/PATH
- Flask Server Proxy: localhost:80/h3lioviz/metadata -> http://localhost:5000/h3lioviz/metadata
- / ->(302) /h3lioviz (not working?)
- /helioviz ->(302) /h3lioviz (not working?)
- /h3lioviz/[FRONTEND-FILENAME]
- /h3lioviz/[FRONTEND-DIR]/[FRONTEND-FILENAME]
- None of the above conditions were met: serve /h3lioviz/index.html

## Logging

`/var/log/apache2/001-pvw_access.log` & `/var/log/apache2/001-pvw_error.log`