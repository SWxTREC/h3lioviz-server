# H3lioViz

This repository hosts the code to run a Paraview server for
the 3D heliospheric output, from codes such as Enlil and Euhforia.

## Local Development and Testing

To test the code locally, you will need the Docker container that has Paraviewweb installed, the test data files, and this cloned git repository. This will require
Docker and the aws-cli programs to run locally.

[Get Docker](https://docs.docker.com/get-docker/)

[Get aws-cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)

NOTE: If you want to run this without aws specify 'base' as your build target. 

## Building and Running Your Own H3lioviz-Server for swxtrec-cdk
### Setup and Build h3lioviz-server

1. Clone the h3lioviz-server repo
2. Download paraview binaries: `curl https://www.paraview.org/paraview-downloads/download.php\?submit\=Download\&version\=v5.10\&type\=binary\&os\=Linux\&downloadFile\=ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz --output ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz`
3. Extract it into  *./docker/binaries*. `tar -xzvf ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz -C ./docker/binaries/ --strip-components=1`
4. Build: `docker build -t h3lioviz .`
5. Sign in: `aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev`
6. Tag the resulting image as a part of the ECR: `docker tag h3lioviz:latest public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev`
7. Push: `docker push public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev`


TODO: Integrate steps

1. Clone the h3lioviz webapp repository: https://bitbucket.lasp.colorado.edu/scm/webapps/h3lioviz.git
2. 


Forcing running EC2 to use new image: 
```bash
export PATH="$PATH:/usr/local/bin"
/docker/docker-launch.sh
```

1. Download the paraview binaries: `curl https://www.paraview.org/paraview-downloads/download.php\?submit\=Download\&version\=v5.10\&type\=binary\&os\=Linux\&downloadFile\=ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz --output ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz`
2. Extract the tarball into ./docker/binaries: `tar -xzvf ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz -C ./docker/binaries/`
3. Extract it into  *./docker/binaries*. You should see *bin*, *lib*, and *share* directories (Note: The tarball is no longer extracted by the Dockerfile, as that inevitably leads to 2 large layers rather than 1).
4. Run `docker build .`

### Pushing Image for Legacy
This assumes that the image has been built using the above steps.
When pushing the image to dev or prod, use `public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:dev` or `public.ecr.aws/swx-trec/pvw-h3lioviz-osmesa:prod` as the ECR image.

1. Log into aws console and create a new public ECR
2. Sign in: `aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin <your-image>`
3. Tag the resulting image as a part of the ECR: `docker tag h3lioviz:latest <your-image>`
4. Push:  `docker push <your-image>`

### Updating a SWx-Trec-cdk Paraview Instance to Use the New Image

1. Navigate to AWS Session Manager in the console and start a session in the paraview instance
2. Edit the .yaml files in /docker/ to point to your image repository
3. Run `export PATH="$PATH:/usr/local/bin/"`
4. Run `/docker/docker-launch.sh` This will update the docker container and relaunch it.
---

### Adding Additional Packages to pvpython

This works by creating a venv and pip installing the requirements.txt file in pvw/requirements.txt. PV_VENV is then exported to tell pvpython where to look for venv. This step is handled by docker/scripts/server.sh. 

1. Add the packages to pvw/requirements.txt
2. For any package that includes `import paraview.simple` or `from paraview import simple` will need to have `import paraview.web.venv` added to it's imports in order to use the packages listed in pvw/requirements.txt

### Getting the Docker container locally
Use your AWS credentials to authorize yourself to the AWS Registry.

```bash
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 135080795405.dkr.ecr.us-east-1.amazonaws.com
```

Pull the container to your local system. (Note: the `osmesa` container is for CPU rendering, if you want to test GPU rendering you'll need the `egl` container)

```bash
docker pull 135080795405.dkr.ecr.us-east-1.amazonaws.com/enlil:pvw-enlil-osmesa
```

### Test data

Two timesteps are included in the test-data directory. To generate the test data for
the docker container you will need to run the processing code to transform the data,
which will expand the data to ~100 MB in size. This can be included in the container
for distribution, so others don't need to deal with any processing.

```bash
python scripts/process_output.py test-data
```

We will store the actual model output files in git, and then the processing files
will be moved to the container. This allows us to update the processing code and
produce new processing routines without blowing up the repository size.

### Running the package locally

Included test cases:

```bash
docker run -p 0.0.0.0:8080:80 -e SERVER_NAME=127.0.0.1:8080 -e PROTOCOL=ws -it public.ecr.aws/enlil/paraview_web_repo:pvw-h3lioviz-osmesa
```

With your own server modifications and data mounted internally:

```bash
docker run -p 0.0.0.0:8080:80 -e SERVER_NAME=127.0.0.1:8080 -e PROTOCOL=ws -v ${PWD}/pvw:/pvw -v ${PWD}/data/pv-ready-data-NEWHASH:/data/pv-ready-data-9a68f9e9 -it public.ecr.aws/enlil/paraview_web_repo:pvw-h3lioviz-osmesa
```

## Building the Dockerfile

This was built following the [guides](https://github.com/Kitware/paraviewweb/tree/master/tools/docker) from Kitware and their repositories. The `dockerfile/` directory contains all of the content that needs to be
copied into the container.

> **_NOTE:_**  You may need to make the scripts executable if they weren't cloned with those permissions.
> `chmod a+x scripts/`

1. Download the Paraview binary for linux and save it into the `docker/binaries/` directory. [https://www.paraview.org/download/](https://www.paraview.org/download/)
    > **_NOTE:_**  If you have a GPU available use the `egl` version, otherwise use the `osemsa`.

2. Build the image locally with an appropriate tag.

    ```bash
    docker build --rm -t public.ecr.aws/enlil/paraview_web_repo:pvw-h3lioviz-osmesa -f docker/Dockerfile .
    ```

    Optionally, deploy the built image to the ECR.
    This requires a public login from the ECR account to
    verify the credentials. You'll want to remove these after
    pushing the image because they can intefere with other
    AWS Docker processes.

    ```bash
    AWS_PROFILE=swx-trec-legacy aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws
    docker push public.ecr.aws/enlil/paraview_web_repo:pvw-h3lioviz-osmesa
    docker logout public.ecr.aws
    ```

3. Run the image setting the proper environment variables and mounting the proper directories.

    ```bash
    docker-compose --env-file docker/.env.local --file docker/docker-compose.yaml up
    ```

    The container requires several input volumes that contain the data directory (that contains
    a `pv-data-3d.nc` input file) mounted at the `/data` location inside the
    container, and the `pvw` directory from this repository mounted
    at the `/pvw` location. To handle that, you can create a `.env.local` environment variable
    file to point to these data locations on your local system. For instance, your local file may
    look like:

    ```bash
    PVW_BACKEND=/home/code/enlil-3d-server/pvw  
    PVW_DATA=/home/data/pv-ready-data
    ```

## Running the server without Docker

To run the server locally on port 1234 the command would be something like the following.

```bash
/path/to/paraview/bin/pvpython pvw/server/pv_server.py --port 1234
```

The frontend can then use wslink to use a websocket to connect to port 1234.

### Installation

The standard Paraview binaries seem to have issues with missing packages. Creating a
new virtual environment and installing the dependencies manually seemed to be the easiest
way to address that. The below commands worked for me to get a Paraview 5.9.1 server running.
The conda pvpython binary is missing `autobahn` and `wslink`, the latter of which isn't
in the conda package manager, so I had to use pip to install it.

```bash
conda create -n paraview python=3.9 paraview autobahn
pip install wslink
```

Now, you can run the server with

```bash
pvpython pvw/server/app_server.py --port 1234 --dir /path/to/<pv-data-3d.nc>
```

where the port is `1234` for local development, and the path to the input file is
given to the `--dir` command line argument.


## Repository Details

### How the container is structured
The container's entrypoint is the server.sh script, which calls start.sh.

 - Because the network adaptor is bridged, the docker container inherits the ec2 instance's permissions
