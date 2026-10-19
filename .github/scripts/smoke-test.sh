#!/usr/bin/env bash
# Prove an image serves a ParaView session before anything publishes it.
#
# Starts the container with only port 80 published: no data, no AWS
# credentials, and none of the EC2's environment (SERVER_NAME, PROTOCOL,
# S3_BUCKET_NAME) or its /data mount, so the session URL rewrite and the run
# metadata routes are not exercised. It waits for the launcher behind Apache,
# then asks it for a session exactly as the frontend does. The launcher answers
# a POST only after the new pvpython process prints its ready line ("Starting
# factory" in pvw/launcher/config.json), so a session URL in the answer means
# ParaView started and pvw/server/app_server.py imported and began serving.
# Last it asks for the metadata server's health route, which the load
# balancer's health check uses, so an image whose flask app does not start
# fails here too.
#
# usage: .github/scripts/smoke-test.sh <image> [host port, default 8080]

set -euo pipefail

IMAGE="$1"
PORT="${2:-8080}"
NAME="h3lioviz-smoke-$$"
# The trailing slash matters: the frontend posts to /h3lioviz/paraview/, and
# without it ParaView 5.10's launcher answers 404
LAUNCHER="http://localhost:${PORT}/h3lioviz/paraview/"

dump_logs() {
    echo "::group::container output"
    docker logs "${NAME}" 2>&1 | tail -100 || true
    echo "::endgroup::"
    echo "::group::launcher logs"
    docker exec "${NAME}" sh -c 'tail -n 60 /data/launcher/log/*' 2>&1 || true
    echo "::endgroup::"
}
cleanup() { docker rm -f "${NAME}" >/dev/null 2>&1 || true; }
trap cleanup EXIT

docker run -d --platform linux/amd64 --name "${NAME}" -p "${PORT}:80" "${IMAGE}" >/dev/null

# A GET is not a session request, so the launcher refuses it: wslink answers
# 400 on ParaView 5.10 and 405 on 6.x. Either one means Apache is proxying to a
# running launcher; anything else (000 while Apache starts, 502/503 while the
# launcher starts) means not yet.
code=000
for _ in $(seq 1 90); do
    code=$(curl -s -o /dev/null -w '%{http_code}' -m 5 "${LAUNCHER}" || true)
    if [ "${code}" = 400 ] || [ "${code}" = 405 ]; then break; fi
    sleep 2
done
if [ "${code}" != 400 ] && [ "${code}" != 405 ]; then
    echo "launcher never came up: last GET answered ${code}"
    dump_logs
    exit 1
fi
echo "launcher is up (GET answered ${code})"

# The launcher waits up to its own 45 s timeout for the session's ready line
answer=$(curl -s -m 90 -w '\n%{http_code}' -X POST -H 'Content-Type: application/json' \
    -d '{"application":"visualizer"}' "${LAUNCHER}" || true)
status=$(printf '%s\n' "${answer}" | tail -n 1)
body=$(printf '%s\n' "${answer}" | sed '$d')
echo "session request answered ${status}: ${body}"
if [ "${status}" != 200 ] || ! printf '%s' "${body}" | grep -q '"sessionURL"'; then
    echo "the launcher did not hand out a session"
    dump_logs
    exit 1
fi
# The flask app starts in the background (docker/scripts/server.sh), so give it
# the same kind of wait as the launcher
health=000
for _ in $(seq 1 30); do
    health=$(curl -s -o /dev/null -w '%{http_code}' -m 5 "http://localhost:${PORT}/h3lioviz/metadata/health" || true)
    if [ "${health}" = 200 ]; then break; fi
    sleep 2
done
if [ "${health}" != 200 ]; then
    echo "the metadata server's health route answered ${health}, not 200"
    dump_logs
    exit 1
fi
echo "metadata health route answered 200"
echo "smoke test passed: ${IMAGE} started a ParaView session and the metadata server"
