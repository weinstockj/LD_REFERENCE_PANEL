#!/bin/bash

# check_requirement.sh
# Script to verify system dependencies for LD Reference Panel Generation pipeline

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local status=$1
    local message=$2
    if [ "$status" = "OK" ]; then
        echo -e "${GREEN}[OK]${NC} $message"
    elif [ "$status" = "WARN" ]; then
        echo -e "${YELLOW}[WARN]${NC} $message"
    elif [ "$status" = "ERROR" ]; then
        echo -e "${RED}[ERROR]${NC} $message"
    fi
}

# Function to check if command exists
check_command() {
    local cmd=$1
    local required_version=$2
    local current_version=""
    
    if command -v "$cmd" >/dev/null 2>&1; then
        case "$cmd" in
            "plink2")
                current_version=$(plink2 --version 2>&1 | head -n1 | grep -oE 'v[0-9]+\.[0-9]+[a-zA-Z0-9.]*' | head -n1)
                ;;
            "bcftools")
                current_version=$(bcftools --version 2>/dev/null | head -n1 | grep -oE '[0-9]+\.[0-9]+[0-9.]*')
                ;;
            "snakemake")
                current_version=$(snakemake --version 2>/dev/null)
                ;;
            "python3")
                current_version=$(python3 --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
                ;;
            *)
                current_version="unknown"
                ;;
        esac
        
        if [ -n "$required_version" ] && [ "$current_version" != "unknown" ]; then
            print_status "OK" "$cmd found (version: $current_version, required: $required_version)"
        else
            print_status "OK" "$cmd found (version: $current_version)"
        fi
        return 0
    else
        print_status "ERROR" "$cmd not found"
        return 1
    fi
}

# Function to check Python packages with uv
check_python_package() {
    local package=$1
    local required_version=$2
    
    # First check if uv is available and try to use it
    if command -v uv >/dev/null 2>&1; then
        if uv run python -c "import $package" >/dev/null 2>&1; then
            local version=$(uv run python -c "import $package; print($package.__version__)" 2>/dev/null || echo "unknown")
            if [ -n "$required_version" ] && [ "$version" != "unknown" ]; then
                print_status "OK" "Python package $package found via uv (version: $version, required: $required_version)"
            else
                print_status "OK" "Python package $package found via uv (version: $version)"
            fi
            return 0
        else
            print_status "ERROR" "Python package $package not found in uv environment"
            return 1
        fi
    else
        # Fallback to system python if uv is not available
        if python3 -c "import $package" >/dev/null 2>&1; then
            local version=$(python3 -c "import $package; print($package.__version__)" 2>/dev/null || echo "unknown")
            if [ -n "$required_version" ] && [ "$version" != "unknown" ]; then
                print_status "OK" "Python package $package found (version: $version, required: $required_version)"
            else
                print_status "OK" "Python package $package found (version: $version)"
            fi
            return 0
        else
            print_status "ERROR" "Python package $package not found"
            return 1
        fi
    fi
}

# Function to check file existence
check_file() {
    local file=$1
    local description=$2
    
    if [ -f "$file" ]; then
        print_status "OK" "$description found: $file"
        return 0
    else
        print_status "ERROR" "$description not found: $file"
        return 1
    fi
}

# Function to check directory existence
check_directory() {
    local dir=$1
    local description=$2
    
    if [ -d "$dir" ]; then
        print_status "OK" "$description found: $dir"
        return 0
    else
        print_status "WARN" "$description not found: $dir"
        return 1
    fi
}

echo "=============================================="
echo "LD Reference Panel - Dependency Check"
echo "=============================================="

# Check system requirements
echo -e "\n1. Checking system requirements..."
check_command "plink2" "v2.00a5.9LM"
check_command "bcftools" "1.15"
check_command "uv"
check_command "python3" ">=3.11"

# Check Python packages
echo -e "\n2. Checking Python packages..."
check_python_package "pandas" ">=2.3.0"
check_python_package "snakemake" ">=9.6.0"
check_python_package "yaml" 

# Check for SLURM (optional but recommended)
echo -e "\n3. Checking job scheduler..."
if command -v "sbatch" >/dev/null 2>&1; then
    print_status "OK" "SLURM found (sbatch available)"
else
    print_status "WARN" "SLURM not found - pipeline can run locally but performance may be limited"
fi

# Check configuration files
echo -e "\n4. Checking configuration files..."
check_file "config.yaml" "Configuration file"
check_file "Snakefile" "Snakefile"

# Check if pyrho data directory exists
echo -e "\n5. Checking data dependencies..."
check_directory "data/pyrho" "Pyrho LD blocks data directory"

# Summary
echo -e "\n=============================================="
echo "Dependency check completed!"
echo "=============================================="

# Check if any critical dependencies are missing
CRITICAL_MISSING=0

if ! command -v plink2 >/dev/null 2>&1; then
    CRITICAL_MISSING=1
fi
if ! command -v bcftools >/dev/null 2>&1; then
    CRITICAL_MISSING=1
fi
# Check Python packages with uv if available, otherwise fallback to system python
if command -v uv >/dev/null 2>&1; then
    if ! uv run python -c "import snakemake" >/dev/null 2>&1; then
        CRITICAL_MISSING=1
    fi
    if ! uv run python -c "import pandas" >/dev/null 2>&1; then
        CRITICAL_MISSING=1
    fi
else
    if ! python3 -c "import snakemake" >/dev/null 2>&1; then
        CRITICAL_MISSING=1
    fi
    if ! python3 -c "import pandas" >/dev/null 2>&1; then
        CRITICAL_MISSING=1
    fi
fi

if [ $CRITICAL_MISSING -eq 1 ]; then
    echo -e "\n${RED}CRITICAL DEPENDENCIES MISSING!${NC}"
    echo "Please install missing dependencies before running the pipeline."
    exit 1
else
    echo -e "\n${GREEN}All critical dependencies satisfied!${NC}"
    echo "Pipeline should run successfully."
    exit 0
fi