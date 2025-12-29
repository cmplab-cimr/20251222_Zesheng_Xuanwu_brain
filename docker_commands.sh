#!/bin/bash

## Function to display usage information
usage() {
  cat << EOF
Usage: bash $0 [OPTIONS] [ARGUMENTS] COMMANDS

A demo script showing advanced getopt features.

Commands:
  start       Start docker
  stop        Stop docker
  log         View docker logs
  result      Execute result process
  build       Build image

Options:
  -p, --parallel NUM     Parallel number for normal slsqp solver. Invalid for tf2 solver (default: 20)
  -s, --suffix STR       Suffix for running (default: none)
  -i, --dataset_id STR   Dataset ID label for running (default: none)
  -h, --help             Display this help message

Examples:
  bash $0 start
  bash $0 -s 1 -i abc start
  bash $0 stop
EOF
  exit 1
}


## Function to start container for flux_analysis and result_process (when parallel_param == 'result')
start_docker() {
  local image_name="$1"
  local project_folder="$2"
  local container_name="$3"
  local code_parent_direct="$4"
  local parallel_param="$5"

  local code_direct="$code_parent_direct/$project_folder"
  echo "Remove old container $container_name and start new in current direct..."
  echo "Code direct: $code_direct"
  docker container rm $container_name

  if [ "$parallel_param" = "result"  ]; then
    echo "Running result_process..."
  else
    echo "Running flux_analysis..."
  fi
  docker run -itd --shm-size=1G --name $container_name -v $code_parent_direct:$code_parent_direct -w $code_direct \
      -e SUFFIX=$SUFFIX -e DATASET_ID=$DATASET_ID \
      $image_name:latest python_main_execute.sh $parallel_param
  EXIT_CODE=$?

  exit $EXIT_CODE
}

## Function to stop container
stop_docker() {
  local container_name="$1"
  echo "Stop container $container_name..."
  docker container stop $container_name
}

## Function to view container log
docker_log() {
  local container_name="$1"
  docker logs -t $container_name | vi -
}

## Function to build image
build_image() {
  local image_name="$1"
  echo "Build image $image_name..."
  docker build . -t $image_name
}

## Parse command-line options
OPTS=$(getopt -o p:s:i:h --long parallel:,suffix:,dataset_id:,help -- "$@")

if [ $? -ne 0 ]; then
  echo "Failed to parse options"
  usage
fi

## Reset the positional parameters to the parsed options
eval set -- "$OPTS"

## Initialize variables with default values
CURRENT_PATH=$PWD
PROJECT_FOLDER=`basename $CURRENT_PATH`
CONTAINER_NAME=${PROJECT_FOLDER}__${USER}
COMMON_SCRIPT_NAME=`dirname $CURRENT_PATH`
IMAGE_NAME="multi_organ_mfa"
PARALLEL_NUM=20
SUFFIX=""
DATASET_ID=""
MODE=""

## Process the options
while true; do
  case "$1" in
    -p | --parallel)
      PARALLEL_NUM="$2"
      echo "Parallel: $PARALLEL_NUM"
      shift 2
      ;;
    -s | --suffix)
      SUFFIX="$2"
      echo "Suffix: '$SUFFIX'"
      shift 2
      ;;
    -i | --dataset_id)
      DATASET_ID="$2"
      echo "Dataset ID: '$DATASET_ID'"
      shift 2
      ;;
    -h | --help)
      usage
      ;;
    --)
      if [ "$2" = "start" ] || [ "$2" = "stop" ] || [ "$2" = "log" ] || [ "$2" = "result" ] || [ "$2" = "build" ]; then
        MODE="$2"
      else
        echo "Error: Invalid command '$2'. Must be 'start', 'stop', 'log' or 'build'."
        usage
      fi
      break
      ;;
    *)
      echo "Error: Invalid parameter '$1' and '$2'. Internal error!"
      exit 1
      ;;
  esac
done

## Check if required options are provided
if [ -z "$MODE" ]; then
  echo "Error: Command must be specified."
  usage
fi

export SUFFIX="$SUFFIX"
export DATASET_ID="$DATASET_ID"

## Run the main execution
if [ "$MODE" = "start"  ]; then
  start_docker $IMAGE_NAME $PROJECT_FOLDER $CONTAINER_NAME $COMMON_SCRIPT_NAME $PARALLEL_NUM
elif [ "$MODE" = "stop"  ]; then
  stop_docker $CONTAINER_NAME
elif [ "$MODE" = "log"  ]; then
  docker_log $CONTAINER_NAME
elif [ "$MODE" = "result"  ]; then
  start_docker $IMAGE_NAME $PROJECT_FOLDER $CONTAINER_NAME $COMMON_SCRIPT_NAME $MODE
elif [ "$MODE" = "build"  ]; then
  build_image $IMAGE_NAME
else
  echo "Error: Invalid command '$MODE'. Must be 'start', 'stop', 'log' or 'build'."
fi
EXIT_CODE=$?

exit $EXIT_CODE