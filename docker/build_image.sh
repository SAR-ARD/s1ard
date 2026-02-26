#!/bin/bash

### ----------------------------------------- ###
### ----------------- Setup ----------------- ###
### ----------------------------------------- ###

# Read flags from CLI/terminal
while getopts 'sn' OPTION; do
  case "$OPTION" in
    s)
        SNAP="true"
        ;;
    n)
        CACHE_CHECK="--no-cache"
        ;;
    ?)
        echo -e "\nScript usage: $(basename $0) [-s] [-n]\n" >&2
        echo -e "\t[-s] - install ESA SNAP inside the Docker image\n"
        echo -e "\t[-n] - Add '--no-cache' option to 'docker build' command\n"
        exit 1
        ;;
  esac
done

# Check working directory for docker directory (which contents are needed for the build process)
if [[ ! -d docker ]]; then
    echo -e "'docker' directory not found in the current path!"
    exit 2
fi

# Make sure to use the Docker buildkit for building
# Ref: https://docs.docker.com/build/buildkit/
export DOCKER_BUILDKIT=1

# Set up container image version and name
S1ARD_IMAGE_NAME="s1ard"
CURRENT_GIT_BRANCH=`git branch --show-current`
CURRENT_GIT_VERSION_LONG=`git describe --tags`
CURRENT_GIT_VERSION_SHORT_NUM=`awk -F "[-]" '{print $1}' <<< $CURRENT_GIT_VERSION_LONG` # Get e.g. 'v0.1.0'
CURRENT_GIT_VERSION_SHORT_COMMIT=`awk -F "[-]" '{print $3}' <<< $CURRENT_GIT_VERSION_LONG` # Get commit SHA (short!)
S1ARD_VERSION="${CURRENT_GIT_VERSION_SHORT_NUM}"


### ------------------------------------------------- ###
### ----------------- Build process ----------------- ###
### ------------------------------------------------- ###

if [[ -z $SNAP ]]; then

    ### ------------------------------------------- ###
    ### ------ Building s1ard-base container ------ ###
    ### ------------------------------------------- ###
    echo -e "\nBuilding s1ard-base container image without SNAP installed...\n"

    if docker build $CACHE_CHECK \
        --target s1ard-base \
        -t ${S1ARD_IMAGE_NAME}:${S1ARD_VERSION} \
        -f Dockerfile . \
    ;then 
        echo "Finished building s1ard-base container!"
    else
        echo "s1ard-base container build failed!"
        exit 10
    fi

else

    ### ------------------------------------------- ###
    ### ------ Building s1ard-snap container ------ ###
    ### ------------------------------------------- ###
    echo -e "\nBuilding s1ard-snap container image with SNAP installed...\n"

    if docker build $CACHE_CHECK \
        --target s1ard-snap \
        -t ${S1ARD_IMAGE_NAME}:${S1ARD_VERSION} \
        -f Dockerfile . \
    ;then 
        echo "Finished building s1ard-snap container!"
    else
        echo "s1ard-snap container build failed!"
        exit 11
    fi

fi