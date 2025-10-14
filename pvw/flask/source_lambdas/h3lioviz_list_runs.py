import json
import boto3
import os
from decimal import *


def lambda_handler(event, context):
    # Get dynamodb table name from environment variable
    table_name = os.environ["TABLE_NAME"]

    # Define dynamodb table
    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table(table_name)

    # Get all items from table
    data = table.scan()["Items"]

    # Convert all Decimal types in table to floats for json conversion
    for item in data:
        for key in item:
            if isinstance(item[key], Decimal):
                item[key] = float(item[key])

    # Return the contents from metadata.json for all runs
    return {
        "statusCode": 200,
        "body": json.dumps(data),
        "headers": {
            "Content-Type": "application/json",
        },
    }
