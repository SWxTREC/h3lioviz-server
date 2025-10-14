import boto3
import os

"""
This function is designed to be invoked by an API Gateway endpoint
The function takes three different routes appended to the API GW URL
/ec2start - To start an EC2
/ec2stop - To stop an EC2
/ec2status - To return both the state and the status of an EC2

Currently, these environmental variables are currently set within the CDK Lambda stack

NOTE: We can use instance tags if we want to start/stop
instances without knowing the ID. We can use a CLI script to
tag instances with their purposes and use this Lambda to start/stop
those instance with tag information in the API request.

"""

# Take in region from the 'input_region' environmental variable
region = os.environ["instance_region"]
# Take in instance id from the 'instance_id' environmental variable
instance = [os.environ["instance_id"]]
ec2 = boto3.client("ec2", region_name=region)


def lambda_handler(event, context):

    # Uncomment to debug function
    print(event)
    print(context)

    # Set the path variable using the path key from the event object
    # An example HTTP API event.json file can be viewed below:
    # https://github.com/awsdocs/aws-lambda-developer-guide/blob/main/sample-apps/nodejs-apig/event-v2.json
    if "rawPath" in event:
        path = event["rawPath"].lower()
    else:
        return {
            "statusCode": 200,
            "body": "No URL path specified",
            "headers": {
                "Content-Type": "text/html",
            },
        }

    # Start by obtaining the instance state to
    ec2_response = ec2.describe_instance_status(
        InstanceIds=instance, IncludeAllInstances=True
    )
    instance_state = ec2_response["InstanceStatuses"][0]["InstanceState"]["Name"]

    # If route path contains 'Start', and instance is stopped,
    # start the EC2 instance and report that the instance is starting
    if "start" in path:
        if instance_state == "stopped":
            ec2.start_instances(InstanceIds=instance)
            content = f"Starting instance {instance}"
        else:
            content = f"{instance} {instance_state}"

    # If route path contains 'Stop', stop the EC2 instance and report that the instance is stopping
    elif "stop" in path:
        ec2.stop_instances(InstanceIds=instance)
        content = f"Stopping instance {instance}"

    # If route path contains 'status', gather the status of the EC2 instance
    elif "status" in path:
        instance_status = ec2_response["InstanceStatuses"][0]["InstanceStatus"][
            "Status"
        ]
        # If the instance is running, it may still be initializing, opened, or have some other status
        # Set the report content to contain this status
        if "running" in instance_state:
            content = f"{instance} {instance_state}, {instance_status}"
        # Otherwise, only report the state of the instance
        else:
            content = f"{instance} {instance_state}"

    # If the path doesn't match any of the above routes, send debug info
    else:
        content = "Invalid option provided via path" + path

    response = {
        "statusCode": 200,
        "body": content,
        "headers": {
            "Content-Type": "text/html",
        },
    }
    return response
