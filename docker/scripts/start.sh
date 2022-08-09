#!/usr/bin/env bash

#
# Patches launcher configuration session url, as well as perhaps any
# additional arguments to pvpython, then restarts the apache webserver
# and starts the launcher in the foreground.  You can optionally pass a
# custom session root url (e.g: "wss://www.example.com") which will be
# used instead of the default.
#
# You can also pass extra arguments after the session url that will be
# provided as extra arguments to pvpython.  In this case, you must also
# pass the session url argument first.
#
# Examples
#
# To just accept the defaults of "ws://localhost":
#
#     ./start.sh
#
# To choose 'wss' and 'www.customhost.com', set the following two environment
# variables before invoking this script:
#
#    export SERVER_NAME="'www.customhost.com"
#    export PROTOCOL="wss"
#
# You can also pass any extra args to pvpython in an environment variable
# as follows:
#
#    export EXTRA_PVPYTHON_ARGS="-dr,--mesa-swr"
#
# Note that the extra args to be passed to pvpython should be separated by
# commas, and no extra spaces are used.
#
# When this script is used as the entrypoint in a Docker image, the environment
# variables can be provided using as many "-e" arguments to the "docker run..."
# command as necessary:
#
#     docker run --rm \
#         -e SERVER_NAME="www.customhost.com" \
#         -e PROTOCOL="wss"
#         -e EXTRA_PVPYTHON_ARGS="-dr,--mesa-swr" \
#         ...
#

ROOT_URL="ws://localhost"
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

# Make sure the apache webserver is running
echo "Starting/Restarting Apache webserver"
service apache2 restart

# Run the pvw launcher in the foreground so this script doesn't end
echo "Starting the wslink launcher"
/opt/paraview/bin/pvpython -m wslink.launcher ${LAUNCHER_PATH}
# python3 -m wslink.launcher ${LAUNCHER_PATH}
