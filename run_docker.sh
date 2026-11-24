#!/bin/bash

# Default values
OS_TYPE="rocky"
FRESH_BUILD=false

# Print script usage instructions
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --ubuntu          Use Ubuntu 24.04 base image"
    echo "  --rocky           Use Rocky Linux 9 base image (Default)"
    echo "  --fresh           Build the image without using cache"
    echo "  -h, --help        Show this help message"
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ubuntu)
            OS_TYPE="ubuntu"
            shift
            ;;
        --rocky)
            OS_TYPE="rocky"
            shift
            ;;
        --fresh)
            FRESH_BUILD=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Set build configuration based on chosen OS
MOUNT_DIR="/EcoSIM"
if [[ "$OS_TYPE" == "ubuntu" ]]; then
    IMAGE_NAME="ubuntu-compiler"
    DOCKERFILE_PATH="docker/ubuntu-compiler.dockerfile"
else
    IMAGE_NAME="rocky-compiler"
    DOCKERFILE_PATH="docker/rocky-compiler.dockerfile"
fi

# Set cache flag
if [ "$FRESH_BUILD" = true ]; then
    CACHE_FLAG="--no-cache"
else
    CACHE_FLAG=""
fi

echo "--- 1. Building the $OS_TYPE Docker image ($IMAGE_NAME) ---"

# Execute the Docker build process
docker build $CACHE_FLAG -f "$DOCKERFILE_PATH" -t "$IMAGE_NAME" .

# Check if build was successful before launching
if [ $? -eq 0 ]; then
    echo "--- 2. Starting the $OS_TYPE container ---"
    docker run -it \
        -v "$(pwd):$MOUNT_DIR" \
        -w "$MOUNT_DIR" \
        "$IMAGE_NAME"
else
    echo "Error: Docker build failed. Container will not start."
    exit 1
fi
