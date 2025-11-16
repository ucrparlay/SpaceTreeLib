#!/bin/bash
# Helper script for running PSI artifact evaluation in Docker

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default values
IMAGE_NAME="ghcr.io/ucrparlay/spacetreelib:latest"
CONTAINER_NAME="psi-artifact-evaluation"
NODE_SIZE="1000000000"

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
usage() {
    cat <<EOF
Usage: $0 [COMMAND] [OPTIONS]

Commands:
    build       Build the Docker image
    pull        Pull the Docker image from GitHub Container Registry
    run         Run the container interactively
    full        Run full artifact evaluation (NODE_SIZE=1000000000)
    shell       Start a bash shell in the container
    stop        Stop the running container
    clean       Remove container and image
    help        Show this help message

Options:
    --data-path PATH    Set the HOST path where data will be stored (required for data generation)
    --node-size SIZE    Set the node size (default: 1000000000)
    --cpus NUM          Limit CPU cores (default: all available)
    --memory SIZE       Memory limit (e.g., 16g, default: no limit)

Examples:
    $0 build
    $0 pull
    $0 run --data-path /mnt/large-disk/experiments
    $0 run --data-path /home/user/data --node-size 100000
    $0 full --data-path /mnt/ssd/experiments --node-size 1000000000 --memory 512g

Note: --data-path specifies where on the HOST machine the generated data will be stored.
      This directory must have sufficient space (100GB+ for full experiments).

Documentation:
    See docs/ folder for detailed guides:
    - docs/QUICKREF.md - Quick command reference
    - docs/DOCKER.md - Detailed Docker instructions
    - docs/ARTIFACT_EVALUATION.md - Complete evaluation guide

EOF
}

# Function to check if Docker is installed
check_docker() {
    if ! command -v docker &>/dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
}

# Function to build the image
build_image() {
    print_info "Building Docker image: $IMAGE_NAME"
    docker build -t "$IMAGE_NAME" .
    print_info "Build complete!"
}

# Function to pull the image
pull_image() {
    print_info "Pulling Docker image from GitHub Container Registry: $IMAGE_NAME"
    docker pull "$IMAGE_NAME"
    print_info "Pull complete!"
}

# Function to create host directories for volume mounts
create_host_dirs() {
    mkdir -p script_ae/data script_ae/logs script_ae/plots
    print_info "Created host directories: script_ae/data, script_ae/logs, script_ae/plots"
}

# Function to run container
run_container() {
    local node_size=$1
    local cpus=$2
    local memory=$3
    local data_path=$4

    if [ -z "$data_path" ]; then
        print_error "Data path is required. Use --data-path to specify where data will be stored on the host."
        print_info "Example: $0 run --data-path /mnt/large-disk/experiments"
        exit 1
    fi

    # Create data directory on host if it doesn't exist
    mkdir -p "$data_path"

    create_host_dirs

    local container_data_path="/data"

    local docker_args=(
        "--cap-add=SYS_NICE"
        # "-u" "$(id -u):$(id -g)"

        "-it"
        "--rm"
        "--name" "$CONTAINER_NAME"
        "-v" "$(pwd)/script_ae/data:/workspace/SpaceTreeLib/script_ae/data"
        "-v" "$(pwd)/script_ae/logs:/workspace/SpaceTreeLib/script_ae/logs"
        "-v" "$(pwd)/script_ae/plots:/workspace/SpaceTreeLib/script_ae/plots"
        "-v" "$data_path:$container_data_path"
        "-e" "DATA_PREFIX=$container_data_path"
        "-e" "NODE_SIZE=$node_size"
    )

    if [ -n "$cpus" ]; then
        docker_args+=("--cpus=$cpus")
        print_info "CPU limit: $cpus cores"
    fi

    if [ -n "$memory" ]; then
        docker_args+=("--memory=$memory")
        print_info "Memory limit: $memory"
    fi

    docker_args+=("$IMAGE_NAME")

    print_info "Host data path: $data_path"
    print_info "Container data path: $container_data_path"
    print_info "NODE_SIZE: $node_size"
    print_info "Live editing enabled - changes on host will reflect in container"
    docker run "${docker_args[@]}"
}

# Function to run full evaluation
run_full_eval() {
    local cpus=$1
    local memory=$2
    local data_path=$3
    local node_size=$4

    if [ -z "$data_path" ]; then
        print_error "Data path is required. Use --data-path to specify where data will be stored on the host."
        print_info "Example: $0 full --data-path /mnt/large-disk/experiments"
        exit 1
    fi

    print_warn "Full evaluation may take several hours!"
    print_info "Running full artifact evaluation with NODE_SIZE=$node_size"

    # Create data directory on host if it doesn't exist
    mkdir -p "$data_path"

    create_host_dirs

    local container_data_path="/data"

    local docker_args=(
        "--cap-add=SYS_NICE"
        "-u" "$(id -u):$(id -g)"

        "-it"
        "--rm"
        "--name" "$CONTAINER_NAME"
        "-v" "$(pwd)/script_ae/data:/workspace/SpaceTreeLib/script_ae/data"
        "-v" "$(pwd)/script_ae/logs:/workspace/SpaceTreeLib/script_ae/logs"
        "-v" "$(pwd)/script_ae/plots:/workspace/SpaceTreeLib/script_ae/plots"
        "-v" "$data_path:$container_data_path"
        "-e" "DATA_PREFIX=$container_data_path"
        "-e" "NODE_SIZE=$node_size"
    )

    if [ -n "$cpus" ]; then
        docker_args+=("--cpus=$cpus")
    fi

    if [ -n "$memory" ]; then
        docker_args+=("--memory=$memory")
    fi

    docker_args+=("$IMAGE_NAME" "bash" "-c" "cd /workspace/SpaceTreeLib/script_ae && ./run.sh $container_data_path $node_size")

    print_info "Host data path: $data_path"
    print_info "Live editing enabled - changes on host will reflect in container"
    docker run "${docker_args[@]}"
    print_info "Full evaluation complete! Results in: results/, logs/, plots/"
}

# Function to open shell in running container
open_shell() {
    if ! docker ps | grep -q "$CONTAINER_NAME"; then
        print_error "Container is not running. Start it first with: $0 run"
        exit 1
    fi

    print_info "Opening shell in running container..."
    docker exec -it "$CONTAINER_NAME" /bin/bash
}

# Function to stop container
stop_container() {
    if docker ps | grep -q "$CONTAINER_NAME"; then
        print_info "Stopping container: $CONTAINER_NAME"
        docker stop "$CONTAINER_NAME"
    else
        print_warn "Container is not running"
    fi
}

# Function to clean up
cleanup() {
    print_info "Cleaning up Docker resources..."

    # Stop container if running
    if docker ps | grep -q "$CONTAINER_NAME"; then
        docker stop "$CONTAINER_NAME"
    fi

    # Remove stopped containers
    docker container prune -f

    # Remove image
    if docker images | grep -q "$IMAGE_NAME"; then
        docker rmi "$IMAGE_NAME"
        print_info "Removed image: $IMAGE_NAME"
    fi

    print_info "Cleanup complete!"
}

# Main script
check_docker

# Parse arguments
COMMAND="${1:-help}"
shift || true

CPUS=""
MEMORY=""
CUSTOM_NODE_SIZE=""
CUSTOM_DATA_PATH=""

while [[ $# -gt 0 ]]; do
    case $1 in
    --data-path)
        CUSTOM_DATA_PATH="$2"
        shift 2
        ;;
    --node-size)
        CUSTOM_NODE_SIZE="$2"
        shift 2
        ;;
    --cpus)
        CPUS="$2"
        shift 2
        ;;
    --memory)
        MEMORY="$2"
        shift 2
        ;;
    *)
        print_error "Unknown option: $1"
        usage
        exit 1
        ;;
    esac
done

# Execute command
case $COMMAND in
build)
    build_image
    ;;
pull)
    pull_image
    ;;
run)
    NODE_SIZE_TO_USE="${CUSTOM_NODE_SIZE:-$NODE_SIZE}"
    run_container "$NODE_SIZE_TO_USE" "$CPUS" "$MEMORY" "$CUSTOM_DATA_PATH"
    ;;
full)
    NODE_SIZE_TO_USE="${CUSTOM_NODE_SIZE:-$NODE_SIZE}"
    run_full_eval "$CPUS" "$MEMORY" "$CUSTOM_DATA_PATH" "$NODE_SIZE_TO_USE"
    ;;
shell)
    open_shell
    ;;
stop)
    stop_container
    ;;
clean)
    cleanup
    ;;
help | --help | -h)
    usage
    ;;
*)
    print_error "Unknown command: $COMMAND"
    usage
    exit 1
    ;;
esac
