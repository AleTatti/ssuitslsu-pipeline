#!/usr/bin/env bash

set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Setup Script for SSU/ITS/LSU Pipeline Environments
# ──────────────────────────────────────────────────────────────────────────────

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local status=$1
    local message=$2
    case $status in
        "info")  echo -e "${GREEN}[INFO]${NC} $message" ;;
        "warn")  echo -e "${YELLOW}[WARN]${NC} $message" ;;
        "error") echo -e "${RED}[ERROR]${NC} $message" ;;
    esac
}

# Check if conda is available
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_status "error" "Conda is not installed or not in PATH"
        exit 1
    fi
    print_status "info" "Found conda: $(conda --version)"
}

# Create conda environment from yaml file
create_env() {
    local yaml_file=$1
    local env_name=$(grep "^name:" "$yaml_file" | awk '{print $2}')

    print_status "info" "Checking environment: $env_name"

    # Check if environment already exists
    if conda env list | grep -q "^${env_name} \|^${env_name}$"; then
        print_status "warn" "Environment '$env_name' already exists"
        read -p "Do you want to update it? (y/N) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_status "info" "Updating environment '$env_name'..."
            conda env update -f "$yaml_file" --prune
        else
            print_status "info" "Skipping '$env_name'"
        fi
    else
        print_status "info" "Creating environment '$env_name'..."
        conda env create -f "$yaml_file"
    fi
}

# Main setup function
main() {
    print_status "info" "Starting SSU/ITS/LSU Pipeline Environment Setup"
    echo

    # Check conda installation
    check_conda
    echo

    # Find the envs directory
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
    ENVS_DIR="$PROJECT_ROOT/envs"

    if [[ ! -d "$ENVS_DIR" ]]; then
        print_status "error" "Environment directory not found: $ENVS_DIR"
        exit 1
    fi

    # List all environment files
    print_status "info" "Found environment files in: $ENVS_DIR"
    ls -1 "$ENVS_DIR"/*.yaml 2>/dev/null | while read -r yaml_file; do
        echo "  - $(basename "$yaml_file")"
    done
    echo

    # Create environments
    ENV_FILES=(
        "$ENVS_DIR/utils.yaml"
        "$ENVS_DIR/fastp.yaml"
        "$ENVS_DIR/itsx.yaml"
        "$ENVS_DIR/chimera.yaml"
        "$ENVS_DIR/taxo.yaml"
        "$ENVS_DIR/mafft.yaml"
        "$ENVS_DIR/iqtree.yaml"
        "$ENVS_DIR/mapping.yaml"
        "$ENVS_DIR/megahit.yaml"
        "$ENVS_DIR/spades.yaml"
    )

    for yaml_file in "${ENV_FILES[@]}"; do
        if [[ -f "$yaml_file" ]]; then
            create_env "$yaml_file"
            echo
        else
            print_status "warn" "File not found: $yaml_file"
        fi
    done

    # Final summary
    echo
    print_status "info" "Setup complete! Created/updated environments:"
    conda env list | grep "ssuitslsu-" || true
    echo
    print_status "info" "You can now run the pipeline with:"
    echo "  ./scripts/run_pipeline.sh [config_file]"
}

# Run the main function
main "$@"
