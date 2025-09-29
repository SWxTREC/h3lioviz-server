#!/usr/bin/env bash

cloud_vars=$(aws cloudformation describe-stacks --stack-name h3lioviz --query "Stacks[0].Outputs" --output json) 
if [ $? -ne 0 ] || grep -q -i error; then
    echo "$cloud_vars"
    echo "Could not retrieve cloud variables"
    exit 1
fi 

cloud_vars_formatted=$(echo "$cloud_vars" | jq 'map({(.ExportName): .OutputValue}) | add')
table_name=$(echo "$cloud_vars_formatted" | jq '.H3liovizRunUsageTableName')
if [ $? -ne 0 ]; then
    echo "$table_name"
    exit 1
fi
insert_output=$(aws dynamodb --table-name "$table_name" --item "$1")
if [ $? -ne 0 ] || grep -q -i error; then
    echo "$insert_output"
    exit 1
fi
