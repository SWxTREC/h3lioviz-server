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
# Can be egl (GPU) or osmesa (CPU)
ARG RENDERER=osmesa

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
        libpciaccess0 \
        python3.12-venv \
        git && \
        rm -rf /var/lib/apt/lists/*

RUN curl -L "https://github.com/peak/s5cmd/releases/download/v2.1.0/s5cmd_2.1.0_Linux-64bit.tar.gz" -o /tmp/s5cmd.tar.gz && \ 
    tar -xvzf /tmp/s5cmd.tar.gz -C /usr/local/bin && \
    rm -rf /tmp/*

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf ./aws awscliv2.zip 

# NOTE: We are using the osmesa build here.  If you want to use EGL (GPU) rendering,
# you'll need to change the URL to point to the EGL build and ensure that your
# base image has the necessary EGL libraries installed.
# For EGL builds, consider using the nvidia/opengl base images.
RUN mkdir -p /opt/paraview && \
    curl -fSL "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v6.1&type=binary&os=Linux&downloadFile=ParaView-6.1.1-MPI-Linux-Python3.12-x86_64.tar.gz" -o /tmp/paraview.tar.gz && \
    tar -xzf /tmp/paraview.tar.gz -C /opt/paraview --strip-components=1 && \
    rm -f /tmp/paraview.tar.gz

# Separate pip build to prevent having to re-install dependencies on every pvw code change.
COPY pvw/requirements.txt /pvw/
RUN python3 -m venv /pvw/venv && source /pvw/venv/bin/activate && pip3 install -r /pvw/requirements.txt --upgrade && deactivate

# The venv is necessary to access pip in this way 
RUN python3 -m venv /pvw/server/venv && source /pvw/server/venv/bin/activate && python3 -m pip install --target /pvw/server/wslink-dependencies wslink

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