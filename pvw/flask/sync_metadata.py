"""Synchronize DynamoDB run metadata with S3 contents."""
import json
import os
from decimal import Decimal
import boto3
import logging
import re

from helpers import _get_env_var

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def get_runs_in_dynamo(ddb_client, table_name): 
    """Returns a set of run ids currently in Dynamo DB"""
    run_ids = set() # Using a set for efficient membership testing
    paginator = ddb_client.get_paginator('scan')
   
    logger.info(f"Getting run ids from DynamoDB table {table_name}")
    try:
        for page in paginator.paginate(
            TableName=table_name,
            Select='SPECIFIC_ATTRIBUTES',
            ProjectionExpression='run_id',
        ):
            # Extract run_ids from the page
            for item in page['Items']:
                if 'run_id' in item:
                    run_ids.add(item['run_id']['S'])
    except Exception as e:
        logger.error(f"Error, failed to scan DynamoDB table {table_name}: {e}")
        raise e   

    return run_ids

def get_runs_in_s3(s3_client, bucket_name):
    """Returns a set of run ids currently in the 'data/h3lioviz/' prefix of the given S3 bucket"""
    run_ids = set()
    paginator = s3_client.get_paginator("list_objects_v2")

    logger.info(f"Getting run ids from S3 bucket {bucket_name}")
    try:
        for page in paginator.paginate(Bucket=bucket_name, Prefix="data/h3lioviz"):
            for obj in filter(lambda x: "Key" in x and x["Key"].endswith("metadata.json"), page.get("Contents", [])):

                # Extract run id from object key
                pattern = r"data/h3lioviz/pv-ready-data-(?P<run_id>\d[^/]*)/metadata.json"
                match = re.search(pattern, obj["Key"])

                if match:
                    run_id = match.group('run_id')
                    run_ids.add(run_id)
                else:
                    logger.warning(f"Pattern did not match for metadata file {obj['Key']}, unable to extract run_id")
                    continue 

    except Exception as e:
        logger.error(f"Error, failed to list object from s3 bucket {bucket_name}: {e}")
        raise e

    return run_ids


def sync_metadata():
    # Get S3 bucket name and DynamoDB table name from environment variables
    bucket_name = _get_env_var("S3_BUCKET_NAME")
    ddb_table_name = _get_env_var("TABLE_NAME")
    
    # Create objects to communicate with S3 and dynamodb
    s3 = boto3.client("s3")
    ddb_client = boto3.client('dynamodb')
    dynamodb = boto3.resource("dynamodb")
    runs_table = dynamodb.Table(ddb_table_name)

    logger.info(f"Syncing bucket {bucket_name}, and table {ddb_table_name}")

    # Get set of runs currently in dynamo db
    try:
        ids_in_ddb = get_runs_in_dynamo(ddb_client, ddb_table_name)
    except Exception as e:
        logger.error(f"Error, failed to get runs in DynamoDB: {e}")
        raise e
    
    # Get set of runs currently in s3
    try:
        ids_in_s3 = get_runs_in_s3(s3, bucket_name)
    except Exception as e:
        logger.error(f"Error, failed to get runs in s3: {e}")
        raise e

    s3_only_runs = ids_in_s3 - ids_in_ddb
    ddb_only_runs = ids_in_ddb - ids_in_s3

    # Delete runs from DynamoDB that are no longer in s3
    for id in ddb_only_runs:
        try:
            response = runs_table.delete_item(
                Key={"run_id": id}
            )
            logger.info(f"Deleted run {id} from {ddb_table_name}")
        except Exception as e:
            logger.error(f"Error deleting run {id} from {ddb_table_name}: {e}")
            continue 

    # Add runs from s3 that are missing in DynamoDB
    for id in s3_only_runs:
        # extract all that data into a dict and send to ddb
        metadata_path=f"data/h3lioviz/pv-ready-data-{id}/metadata.json"

        # Read metadata.json file into python dictionary
        try:
            data = json.loads(
                s3.get_object(Bucket=bucket_name, Key=metadata_path)["Body"].read(),
                parse_float=Decimal,
            )
        except Exception as e:
            logger.error(f"Error, failed to parse object {metadata_path}: {e}")
            continue 

        data["last_accessed"] = 0 # Add last_accessed value and set to 0
        logger.info(f"Adding metadata for run {id} into table {ddb_table_name}")
        try:
            runs_table.put_item(Item=data)
        except Exception as e:
            logger.error(f"Error, failed to insert metadata for run {id} into {ddb_table_name}: {e}")