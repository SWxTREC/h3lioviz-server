ARG BASE_IMAGE=ubuntu:20.04
# Setting the below argument to false disables installing AWS packages
ARG AWS_ENABLED=false

# ARG BASE_IMAGE=nvidia/opengl:1.0-glvnd-devel-ubuntu20.04
FROM ${BASE_IMAGE} AS base

# Install using bash rather than sh. Allows for sourcing the pip venv. 
SHELL ["/bin/bash", "-c"]

USER root

# Need to force noninteractive for apt-get updates
ARG DEBIAN_FRONTEND=noninteractive
# Can be egl (GPU) or osmesa (CPU)
ARG RENDERER=osmesa
ARG SESSION_MANAGER_URL="https://paraview-web.noaa-demo.swx-trec.com/h3lioviz/paraview/"
ARG API_URL="https://paraview-web.noaa-demo.swx-trec.com/h3lioviz/metadata/"



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
        python3.9 \
        python3.9-venv \
        git \
        nodejs \
        npm && \
        rm -rf /var/lib/apt/lists/*

# Instructions for installing node for  https://deb.nodesource.com
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash -
RUN sudo apt-get install -y nodejs

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
    curl -fSL "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.10&type=binary&os=Linux&downloadFile=ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz" -o /tmp/paraview.tar.gz && \
    tar -xzf /tmp/paraview.tar.gz -C /opt/paraview --strip-components=1 && \
    rm -f /tmp/paraview.tar.gz

# Separate pip build to prevent having to re-install dependencies on every pvw code change.
COPY pvw/requirements.txt /pvw/
RUN python3.9 -m venv /pvw/venv && source /pvw/venv/bin/activate && pip3 install -r /pvw/requirements.txt --upgrade && deactivate

RUN groupadd proxy-mapping && \
    groupadd pvw-user && \
    useradd --system -g pvw-user -G proxy-mapping -s /sbin/nologin pvw-user && \
    usermod -a -G proxy-mapping www-data && \
    useradd admin && echo "admin:admin" | chpasswd && adduser admin sudo && \
    mkdir -p /opt/launcher/log && \
    chown -R pvw-user:pvw-user /opt/launcher && \
    mkdir -p /opt/paraviewweb/scripts && \
    touch /opt/launcher/proxy-mapping.txt && \
    chown pvw-user:proxy-mapping /opt/launcher/proxy-mapping.txt && \
    chmod 660 /opt/launcher/proxy-mapping.txt

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

RUN mkdir /data

# NOTE: Having the frontend build 
RUN if [ ! -d "/pvw/www/h3lioviz" ]; then \
    git clone https://github.com/SWxTREC/h3lioviz.git /h3lioviz && \
    cd /h3lioviz && \
    sed -i "s|[[:space:]]*sessionManagerURL:.*|sessionManagerURL: \"${SESSION_MANAGER_URL}\",|" src/environments/environment.dev.ts &&\
    sed -i "s|[[:space:]]*api:.*|api: \"${API_URL}\",|" src/environments/environment.dev.ts && \
    npm install && \
    npm run build:dev && \
    cp -r dist/h3lioviz /pvw/www/h3lioviz && \
    cd /; \
else \
    echo "pvw/www/h3lioviz found. Skipping frontend h3lioviz build."; \
fi

# Start the container.  If we're not running this container, but rather are
# building other containers based on it, this entry point can/should be
# overridden in the child container.  In that case, use the "start.sh"
# script instead, or you can provide a custom one.
ENTRYPOINT ["/opt/paraviewweb/scripts/server.sh"]