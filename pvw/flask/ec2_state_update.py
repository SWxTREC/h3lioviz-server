"""Helpers for managing the EC2 lifecycle endpoints."""

from __future__ import annotations

import os
from typing import Dict

import boto3
from botocore.exceptions import ClientError


def _get_client():
    region = os.environ.get("instance_region")
    if not region:
        raise RuntimeError("Missing required environment variable instance_region")
    return boto3.client("ec2", region_name=region)


def _get_instance_id() -> str:
    instance_id = os.environ.get("instance_id")
    if not instance_id:
        raise RuntimeError("Missing required environment variable instance_id")
    return instance_id


def _describe_instance(client, instance_id: str) -> Dict:
    response = client.describe_instance_status(
        InstanceIds=[instance_id], IncludeAllInstances=True
    )
    statuses = response.get("InstanceStatuses", [])
    if not statuses:
        raise RuntimeError(f"Instance {instance_id} not found or has no status")
    return statuses[0]


def start_instance() -> Dict[str, str]:
    client = _get_client()
    instance_id = _get_instance_id()
    status = _describe_instance(client, instance_id)
    current_state = status["InstanceState"]["Name"]

    message: str
    if current_state == "stopped":
        client.start_instances(InstanceIds=[instance_id])
        message = f"Starting instance {instance_id}"
    else:
        message = f"{instance_id} {current_state}"

    return {"message": message, "state": current_state}


def stop_instance() -> Dict[str, str]:
    client = _get_client()
    instance_id = _get_instance_id()
    client.stop_instances(InstanceIds=[instance_id])
    return {"message": f"Stopping instance {instance_id}", "state": "stopping"}


def get_instance_status() -> Dict[str, str]:
    client = _get_client()
    instance_id = _get_instance_id()
    status = _describe_instance(client, instance_id)
    state = status["InstanceState"]["Name"]
    if state == "running":
        instance_status = status["InstanceStatus"]["Status"]
        message = f"{instance_id} {state}, {instance_status}"
    else:
        message = f"{instance_id} {state}"
    return {"message": message, "state": state}
