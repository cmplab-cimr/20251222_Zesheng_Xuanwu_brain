#!/usr/bin/env bash

suffix_input=$SUFFIX
if [ ! $suffix_input ]; then
  suffix_parameter=""
else
  suffix_parameter="-s $suffix_input"
fi

dataset_id_input=$DATASET_ID
if [ ! $dataset_id_input ]; then
  dataset_id_parameter=""
else
  dataset_id_parameter="-i $dataset_id_input"
fi

result=$1
if [ "$result" = "result"  ]; then
  python main.py computation experimental_mfa result_process
else
  while true
  do
    python main.py computation experimental_mfa flux_analysis $suffix_parameter $dataset_id_parameter -p $1
    if [ $? != "3" ]; then
        break
    fi
  done
fi
