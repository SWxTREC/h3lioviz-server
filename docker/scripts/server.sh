#!/usr/bin/env bash

# Note: This is commented out because we don't require packages for paraview.
# Check for a requirements.txt in the mounted root, "/pvw/requirements.txt".
# If we find it, we can pip install those.
# if [ -f "/pvw/requirements.txt" ]; then
#   python3.9 -m venv /pvw/venv
#   source /pvw/venv/bin/activate
#   pip3 install -r /pvw/requirements.txt --upgrade
#   deactivate
#   export PV_VENV="/pvw/venv/"
# fi

# Copy the launcher config into the location where the start script expects
# to find it.  The config may or may not have replacement values in it, if it
# does not, the start script will not change it in any way.  Here we expect
# that the user doing the "docker run ..." has set up an external directory
# containing a "launcher/config.json" filepath and mounts that path as "/pvw".
cp /pvw/launcher/config.json /opt/launcher/config-template.json

# This performs replacements on the launcher-template.json copied into place
# above, based on the presence of environment variables passed with "-e" to the
# "docker run ..." command.
/opt/paraviewweb/scripts/start.sh
