
## Building and Running Your Own H3lioviz-Server for swxtrec-cdk
## Setup h3lioviz-server for Building

1. Navigate to https://www.paraview.org/download/ and download *ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz* which will actually download as *ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar*
2. Copy it into  *./docker/binaries*

## Pushing image for legacy

The bryan-test/h3lioviz public ECR i
1. Sign in: `aws ecr-public get
-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/enlil`
2. Build: `docker build -t h3lioviz .`
3. Tag the resulting image as a part of the ECR: `docker tag h3lioviz:latest public.ecr.aws/enlil/bryan-test/h3lioviz:latest`
4. Push:  `docker push public.ecr.aws/enlil/bryan-test/h3lioviz:latest`

## Updating swx-trec-cdk

1. Navigate to the ec2_construct.py within h3lioviz and change the used ECR to point to `public.ecr.aws/enlil/bryan-test/h3lioviz:latest`
2. Destroy the h3lioviz stack in order to cause the EC2 to be re initialiazed
3. Deploy the h3lioviz stack