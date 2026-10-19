# Ubuntu 24.04: 20.04 is out of standard support, and its Python 3.9 (used for
# flask) is past end of life. 24.04's system Python is 3.12.
ARG BASE_IMAGE=ubuntu:24.04
# Setting the below argument to false disables installing AWS packages
ARG AWS_ENABLED=false
ARG BUILD_FRONTEND=true

# ============================================================================
# Frontend builder stage - builds h3lioviz in a Node.js container
# This stage only rebuilds when the h3lioviz repository changes
# and the BUILD_FRONTEND arg is set to true
# ============================================================================
FROM node:22-slim AS frontend-builder

ARG FRONTEND_ENVIRONMENT="dev"
ARG H3LIOVIZ_VERSION=main

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        unzip \
        ca-certificates && \
        rm -rf /var/lib/apt/lists/*


# Download and build h3lioviz
# The H3LIOVIZ_VERSION arg can be used to bust cache when the repo changes
RUN curl -L "https://github.com/SWxTREC/h3lioviz/archive/refs/heads/${H3LIOVIZ_VERSION}.zip" -o /tmp/h3lioviz.zip && \
    unzip /tmp/h3lioviz.zip -d /tmp/h3lioviz && \
    rm -rf /tmp/h3lioviz.zip

WORKDIR /tmp/h3lioviz/h3lioviz-${H3LIOVIZ_VERSION}

RUN npm install --prefer-offline --no-audit --progress=false && \
    npm rebuild esbuild && \
    npm run build:${FRONTEND_ENVIRONMENT}

# ============================================================================
# Main application stage
# ============================================================================
# ARG BASE_IMAGE=nvidia/opengl:1.0-glvnd-devel-ubuntu20.04
FROM ${BASE_IMAGE} AS base

# Install using bash rather than sh. Allows for sourcing the pip venv. 
SHELL ["/bin/bash", "-c"]

USER root

# Need to force noninteractive for apt-get updates
ARG DEBIAN_FRONTEND=noninteractive
ARG FRONTEND_ENVIRONMENT="dev"

RUN apt-get update && apt-get install -y --no-install-recommends \
        apache2-dev \
        apache2 \
        libapr1-dev \
        libglapi-mesa \
        apache2-utils \
        sudo \
        curl \
        ca-certificates \
        unzip \
        python3 \
        python3-venv \
        git && \
        rm -rf /var/lib/apt/lists/*

RUN curl -L "https://github.com/peak/s5cmd/releases/download/v2.1.0/s5cmd_2.1.0_Linux-64bit.tar.gz" -o /tmp/s5cmd.tar.gz && \ 
    tar -xvzf /tmp/s5cmd.tar.gz -C /usr/local/bin && \
    rm -rf /tmp/*

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf ./aws awscliv2.zip 

# Since ParaView 6.0 there is one Linux binary for every renderer (X, EGL and
# a bundled OSMesa). VTK_DEFAULT_OPENGL_WINDOW below picks OSMesa.
ARG PARAVIEW_VERSION=6.1.1
ARG PARAVIEW_PYTHON=3.12
RUN mkdir -p /opt/paraview && \
    curl -fSL "https://www.paraview.org/files/v${PARAVIEW_VERSION%.*}/ParaView-${PARAVIEW_VERSION}-MPI-Linux-Python${PARAVIEW_PYTHON}-x86_64.tar.gz" -o /tmp/paraview.tar.gz && \
    tar -xzf /tmp/paraview.tar.gz -C /opt/paraview --strip-components=1 && \
    rm -f /tmp/paraview.tar.gz

# System libraries the ParaView 6 binaries link against but do not bundle.
# pvpython will not start without libpciaccess, and without the rest of these
# rendering segfaults, even with OSMesa doing the drawing.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpciaccess0 \
        libgl1 \
        libegl1 \
        libx11-6 \
        libxkbcommon0 \
        libdbus-1-3 \
        libfontconfig1 && \
        rm -rf /var/lib/apt/lists/*

# ParaView 6 no longer bundles wslink, the websocket server every h3lioviz
# session runs on. Install it into a venv on the same Python version as
# pvpython, which loads it through `pvpython --venv /opt/pv-venv`. Stay on
# wslink 1.x: 2.x switched the wire format to msgpack and would need the
# frontend's wslink client upgraded at the same time.
# ParaView 6.2 ships uv in /opt/paraview/bin but 6.1 does not, so take uv
# from its official image, pinned to the version ParaView 6.2.0 bundles.
COPY --from=ghcr.io/astral-sh/uv:0.11.29 /uv /usr/local/bin/uv
ARG WSLINK_VERSION=1.12.4
ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python
RUN uv venv --python ${PARAVIEW_PYTHON} /opt/pv-venv && \
    uv pip install --python /opt/pv-venv/bin/python "wslink==${WSLINK_VERSION}"

# Render with OSMesa as 5.10 did. Left alone, ParaView 6 picks EGL, which
# crashed at pvpython exit in this image.
ENV VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow

# Separate pip build to prevent having to re-install dependencies on every pvw code change.
COPY pvw/requirements.txt /pvw/
RUN python3 -m venv /pvw/venv && source /pvw/venv/bin/activate && pip3 install -r /pvw/requirements.txt --upgrade && deactivate

# The admin user is not given sudo
RUN groupadd proxy-mapping && \
    groupadd pvw-user && \
    useradd --system -g pvw-user -G proxy-mapping -s /sbin/nologin pvw-user && \
    usermod -a -G proxy-mapping www-data && \
    useradd admin && echo "admin:admin" | chpasswd && \
    mkdir -p /opt/launcher/log && \
    chown -R pvw-user:pvw-user /opt/launcher && \
    mkdir -p /opt/paraviewweb/scripts && \
    touch /opt/launcher/proxy-mapping.txt && \
    chown pvw-user:proxy-mapping /opt/launcher/proxy-mapping.txt && \
    chmod 660 /opt/launcher/proxy-mapping.txt

RUN mkdir -p /data/launcher/log

# Copy the apache configuration file into place
COPY docker/config/apache/001-pvw.conf /etc/apache2/sites-available/001-pvw.conf

# Copy the script into place
COPY docker/scripts/* /opt/paraviewweb/scripts/

# Configure the apache web server
RUN a2enmod vhost_alias && \
    a2enmod proxy && \
    a2enmod proxy_http && \
    a2enmod proxy_wstunnel && \
    a2enmod rewrite && \
    a2enmod headers && \
    a2dissite 000-default.conf && \
    a2ensite 001-pvw.conf && \
    a2dismod autoindex -f && \
    apachectl configtest || (cat /var/log/apache2/error.log && exit 1)

# Open port 80 to the world outside the container
EXPOSE 80

# Copy our server release into the container as well
# This can be overridden for local testing with `-v ${PWD}/pvw:/pvw`
COPY pvw /pvw

# Copy the built frontend from the builder stage
# This only rebuilds when the frontend-builder stage changes
FROM base AS frontend-true
COPY --from=frontend-builder /tmp/h3lioviz/h3lioviz-main/dist/h3lioviz /pvw/www/h3lioviz

# Do nothing if BUILD_FRONTEND is set to false
FROM base AS frontend-false

FROM frontend-${BUILD_FRONTEND} AS final 

# Start the container.  If we're not running this container, but rather are
# building other containers based on it, this entry point can/should be
# overridden in the child container.  In that case, use the "start.sh"
# script instead, or you can provide a custom one.
ENTRYPOINT ["/opt/paraviewweb/scripts/server.sh"]