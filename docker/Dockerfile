ARG BASE_IMAGE=ubuntu:24.04

# ARG BASE_IMAGE=nvidia/opengl:1.0-glvnd-devel-ubuntu20.04
FROM ${BASE_IMAGE}

USER root

# Need to force noninteractive for apt-get updates
ARG DEBIAN_FRONTEND=noninteractive
# Can be egl (GPU) or osmesa (CPU)
ARG RENDERER=osmesa

RUN apt-get update && apt-get install -y --no-install-recommends \
        apache2-dev \
        apache2 \
        libapr1-dev \
        libglapi-mesa \
        apache2-utils \
        libx11-6 \
        libgl1 \
        libxcursor-dev \
        python3-pip \
        # libosmesa6-dev?
        # TODO: Look into wslink missing module?
        # /opt/paraview/bin/pvpython-real: Error while finding module specification for 'wslink.launcher' (ModuleNotFoundError: No module named 'wslink')
        sudo && \
    rm -rf /var/lib/apt/lists/*

RUN groupadd proxy-mapping && \
    groupadd pvw-user && \
    useradd --system -g pvw-user -G proxy-mapping -s /sbin/nologin pvw-user && \
    usermod -a -G proxy-mapping www-data && \
    useradd admin && echo "admin:admin" | chpasswd && sudo usermod -aG sudo admin && \
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
    a2dismod autoindex -f

# Open port 80 to the world outside the container
EXPOSE 80

# Copy in the tar file, extract it, then rename it /opt/paraview/...
# COPY binaries/ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit.tar.gz /opt
# RUN cd /opt && \
#     tar -xzvf ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit.tar.gz && \
#     mv ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit paraview

# COPY binaries/ParaView-5.9.1-egl-MPI-Linux-Python3.8-64bit.tar.gz /opt
# RUN cd /opt && \
#     tar -xzvf ParaView-5.9.1-egl-MPI-Linux-Python3.8-64bit.tar.gz && \
#     mv ParaView-5.9.1-egl-MPI-Linux-Python3.8-64bit paraview

# COPY docker/binaries/ParaView-5.10.1-${RENDERER}-MPI-Linux-Python3.9-x86_64.tar.gz /opt
# RUN cd /opt && \
#     tar -xzvf ParaView-5.10.1-${RENDERER}-MPI-Linux-Python3.9-x86_64.tar.gz && \
#     mv ParaView-5.10.1-${RENDERER}-MPI-Linux-Python3.9-x86_64 paraview && \
#     rm ParaView-5.10.1-${RENDERER}-MPI-Linux-Python3.9-x86_64.tar.gz

COPY docker/binaries/ParaView-6.0.0-RC2-MPI-Linux-Python3.12-x86_64.tar.gz /opt
RUN cd /opt && \
    tar -xzvf ParaView-6.0.0-RC2-MPI-Linux-Python3.12-x86_64.tar.gz && \
    mv ParaView-6.0.0-RC2-MPI-Linux-Python3.12-x86_64 paraview && \
    rm ParaView-6.0.0-RC2-MPI-Linux-Python3.12-x86_64.tar.gz

# Copy our server release into the container as well
# This can be overridden for local testing with `-v ${PWD}/pvw:/pvw`
COPY pvw /pvw

# COPY the test data into the container for a default launch case
RUN mkdir /data
COPY test-data/launcher /data/launcher
COPY test-data/pv-ready-data-9a68f9e9/ /data/pv-ready-data-9a68f9e9/

# Need to install this to get wslink installed and loadable as a module.
# NOTE: It is not a part of pvpython, does that matter?
RUN python3 -m pip install "vtk[web]" --break-system-packages

# Start the container.  If we're not running this container, but rather are
# building other containers based on it, this entry point can/should be
# overridden in the child container.  In that case, use the "start.sh"
# script instead, or you can provide a custom one.
ENTRYPOINT ["/opt/paraviewweb/scripts/server.sh"]

# TESTING:
# install vtk[web] for wslink: pip install vtk[web] --break-system-packages
