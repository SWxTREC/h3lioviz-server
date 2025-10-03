ARG BASE_IMAGE=ubuntu:20.04

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
        ca-certificates  \
        curl \
        unzip \
        sudo && \
    rm -rf /var/lib/apt/lists/*

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" && \
    unzip /tmp/awscliv2.zip -d /tmp/ && \
    /tmp/aws/install && \
    rm -rf /tmp/*

RUN curl -L "https://github.com/peak/s5cmd/releases/download/v2.1.0/s5cmd_2.1.0_Linux-64bit.tar.gz" -o /tmp/s5cmd.tar.gz && \ 
    tar -xvzf /tmp/s5cmd.tar.gz -C /usr/local/bin && \
    rm -rf /tmp/*

COPY docker/binaries/ /opt/paraview/


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

# Copy our server release into the container as well
# This can be overridden for local testing with `-v ${PWD}/pvw:/pvw`
COPY pvw /pvw

RUN mkdir /data

# Start the container.  If we're not running this container, but rather are
# building other containers based on it, this entry point can/should be
# overridden in the child container.  In that case, use the "start.sh"
# script instead, or you can provide a custom one.
ENTRYPOINT ["/opt/paraviewweb/scripts/server.sh"]
