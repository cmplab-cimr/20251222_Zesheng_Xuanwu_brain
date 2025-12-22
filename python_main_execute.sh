#!/usr/bin/env bash

suffix_input=$2
if [ ! $suffix_input ]; then
  suffix_parameter=""
else
  suffix_parameter="-s $suffix_input"
fi

result=$1
if [ "$result" = "result"  ]; then
  python main.py computation experimental_mfa result_process
else
  while true
  do
    python main.py computation experimental_mfa flux_analysis $suffix_parameter -p $1
    if [ $? != "3" ]; then
        break
    fi
  done
fi
