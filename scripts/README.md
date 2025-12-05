# Generating H3lioviz Compatible Runs Without Lambda
| Instructions | Commands | Example |
|---|---|---|
|Enter the scripts/ directory|`cd scripts` |
|Create a Python3.12 environment to install dependencies (works with 3.9 as well)|`python3.12 -m venv venv`|
|Activate the Python environment|`source /venv/bin/activate`|
|Install needed dependencies|`pip3 install -r requirements.txt`|
|Process the data by calling `process_output.py` with the path to the directory which contains all the data files. **If you are processing SWPC data** the final directory must match the following regex for run_id to be correctly as run_id is identified from the path. Otherwise a hash of the metadata will be used for the run_id. Regex for SWPC `^.*wsa_enlil_\d{5}\.\d*\.dbqs0`. **Note that the downsampling is configured by default** (8x radius, 2x latitude, and 2x longitude downsample) and you do not need to pass any flags. If you would like to create custom downsamples, you can run `python3 process_output.py -h` to find the appropriate flags.|`python3 process_output.py <path_to_nc_files>`|`python3 process_output.py ~/wsa_enlil_57671.77346046.dbqs01/`|
|Verify the processing was successful by looking for <path_to_nc_files>/pv-ready-data-<run_id>||`ls ~/wsa_enlil_57671.77346046.dbqs01/pv-ready-data-57671/`|
|Upload the run data to the appropriate s3 path. The example uses the aws CLI, but quicker copies of large number of files can be achieved with an open source tool called [s5cmd](https://github.com/peak/s5cmd) which parallelizes uploads.|`aws s3 sync <path_to_nc_files>/pv-ready-data-<run_id>/ s3://h3lioviz.<domain>/data/h3lioviz/pv-ready-data-<run_id>`|`aws s3 sync ~/wsa_enlil_57671.77346046.dbqs01/pv-ready-data-57671 s3://h3lioviz.bryandev.swx-trec.com/data/h3lioviz/pv-ready-data-57671`|
# Getting the New Run Populated into Dynamo DB
| Instructions | Commands | Example |
|---|---|---|
|Log into the paraview ec2|||
|Send a GET http request to the route `/h3lioviz/metadata/syncMetadata`. This will add any runs found in S3 that haven't had their metadata populated into the Dynamo DB table and will also remove any entries in Dynamo DB that do not have an associated metadata.json in S3.|`curl localhost/h3lioviz/metadata/syncMetadata`|
|You should get a status 'null' back. Anything else indicates an error.|
|(Optional) Validate that all the runs in S3 that have a metadata.json file have entries in the dynamodb table and vice versa|
|The run should now be available in the ParaView web application upon refreshing the page|