#!/usr/bin/env bash

# Note: This is commented out because we don't require packages for paraview. 
#       This allows paraview to access the packages in the venv created by the 
#       Dockerfile.
# if [ -f "/pvw/requirements.txt" ]; then
#   export PV_VENV="/pvw/venv/"
# fi

# Copy the launcher config into the location where the start script expects
# to find it. The config may or may not have replacement values in it, if it
# does not, the start script will not change it in any way.  Here we expect
# that the user doing the "docker run ..." has set up an external directory
# containing a "launcher/config.json" filepath and mounts that path as "/pvw".
cp /pvw/launcher/config.json /opt/launcher/config-template.json

# Patches additional parameters into the Paraview-Web config file, which 
# are added to the pvpython startup command. By passing the EXTRA_PVPYTHON_ARGS
# environment variable, we can pass flags to pvpython.

# Available environment variables:
#   - SERVER_NAME: the server name to use for the session url 
#       Gets returned from paraview-web to tell the frontend where to connect to the websocket
#   - PROTOCOL: the protocol to use for the session url 
#       ws -> websocket
#       wss -> websocket secure
#   - EXTRA_PVPYTHON_ARGS: extra arguments to pass to pvpython (comma-separated, no extra spaces)
#       Example: "-dr,--mesa-swr"

ROOT_URL="ws://localhost"  # Default root
REPLACEMENT_ARGS=""

LAUNCHER_TEMPLATE_PATH=/opt/launcher/config-template.json
LAUNCHER_PATH=/opt/launcher/config.json

if [[ ! -z "${SERVER_NAME}" ]] && [[ ! -z "${PROTOCOL}" ]]
then
  ROOT_URL="${PROTOCOL}://${SERVER_NAME}"
fi

if [[ ! -z "${EXTRA_PVPYTHON_ARGS}" ]]
then
  IFS=',' read -ra EXTRA_ARGS <<< "${EXTRA_PVPYTHON_ARGS}"
  for arg in "${EXTRA_ARGS[@]}"; do
    REPLACEMENT_ARGS="${REPLACEMENT_ARGS}\"$arg\", "
  done
fi

INPUT=$(<"${LAUNCHER_TEMPLATE_PATH}")
OUTPUT="${INPUT//"SESSION_URL_ROOT"/$ROOT_URL}"
OUTPUT="${OUTPUT//"EXTRA_PVPYTHON_ARGS"/$REPLACEMENT_ARGS}"
echo -e "$OUTPUT" > "${LAUNCHER_PATH}"

echo "Starting flask webserver"
source /pvw/venv/bin/activate
python3.9 /pvw/flask/main.py >> /data/launcher/log/flask.log 2>&1 &

# Make sure the apache webserver is running
echo "Starting/Restarting Apache webserver"
service apache2 restart

# Run the pvw launcher in the foreground so this script doesn't end
echo "Starting the wslink launcher"
/opt/paraview/bin/pvpython -m wslink.launcher ${LAUNCHER_PATH}
# python3 -m wslink.launcher ${LAUNCHER_PATH}
