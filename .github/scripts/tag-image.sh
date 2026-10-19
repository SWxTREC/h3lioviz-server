#!/usr/bin/env bash
# Point a tag at an image already in the public repository, by digest.
#
# Nothing is pulled or rebuilt. The image's manifest is fetched from the public
# registry without credentials, checked against the digest, and put back with
# the new tag (ECR Public PutImage), so the tag names exactly that image.
# The tag stops naming whatever image it named before; that image stays in the
# repository under any other tags it has.
#
# usage: .github/scripts/tag-image.sh <repository> <sha256:digest> <tag>
# needs: ecr-public:PutImage on the repository (AWS credentials in the environment)

set -euo pipefail

REPOSITORY="$1"
DIGEST="$2"
TAG="$3"
ALIAS="${REGISTRY_ALIAS:-swx-trec}"

if ! [[ "${DIGEST}" =~ ^sha256:[0-9a-f]{64}$ ]]; then
    echo "not a sha256 digest: ${DIGEST}"
    exit 1
fi

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

TOKEN=$(curl -sf "https://public.ecr.aws/token/?scope=repository:${ALIAS}/${REPOSITORY}:pull" \
    | python3 -c 'import json, sys; print(json.load(sys.stdin)["token"])')
curl -sf -o "${WORK}/manifest" -D "${WORK}/headers" \
    -H "Authorization: Bearer ${TOKEN}" \
    -H 'Accept: application/vnd.docker.distribution.manifest.v2+json, application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.list.v2+json, application/vnd.oci.image.index.v1+json' \
    "https://public.ecr.aws/v2/${ALIAS}/${REPOSITORY}/manifests/${DIGEST}" \
    || { echo "no image ${DIGEST} in ${ALIAS}/${REPOSITORY}"; exit 1; }

# The digest is the sha256 of these exact bytes; refuse anything else
GOT="sha256:$(sha256sum "${WORK}/manifest" | cut -d' ' -f1)"
if [ "${GOT}" != "${DIGEST}" ]; then
    echo "the registry returned a manifest with digest ${GOT}, not ${DIGEST}"
    exit 1
fi
MEDIA_TYPE=$(tr -d '\r' < "${WORK}/headers" | awk -F': ' 'tolower($1) == "content-type" {print $2}' | tail -n 1)

aws ecr-public put-image --region us-east-1 \
    --repository-name "${REPOSITORY}" \
    --image-tag "${TAG}" \
    --image-digest "${DIGEST}" \
    --image-manifest-media-type "${MEDIA_TYPE}" \
    --image-manifest "file://${WORK}/manifest" \
    --query 'image.imageId' --output text
echo "${TAG} -> ${DIGEST}"
