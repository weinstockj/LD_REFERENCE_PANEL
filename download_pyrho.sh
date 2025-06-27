#!/bin/bash

# download_pyrho.sh
# Script to download pyrho LD block files from jmacdon/LDblocks_GRCh38 repository

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

# GitHub repository details
REPO_URL="https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data"
DATA_DIR="data/pyrho"

# List of pyrho files to download
PYRHO_FILES=(
    "pyrho_AFR_LD_blocks.bed"
    "pyrho_EAS_LD_blocks.bed"
    "pyrho_EUR_LD_blocks.bed"
    "pyrho_SAS_LD_blocks.bed"
)

echo "=============================================="
echo "Pyrho LD Blocks Data Download"
echo "=============================================="

# Create data directory if it doesn't exist
if [ ! -d "$DATA_DIR" ]; then
    print_status "OK" "Creating data directory: $DATA_DIR"
    mkdir -p "$DATA_DIR"
else
    print_status "OK" "Data directory exists: $DATA_DIR"
fi

# Check if wget or curl is available
DOWNLOAD_CMD=""
if command -v wget >/dev/null 2>&1; then
    DOWNLOAD_CMD="wget"
    print_status "OK" "Using wget for downloads"
elif command -v curl >/dev/null 2>&1; then
    DOWNLOAD_CMD="curl"
    print_status "OK" "Using curl for downloads"
else
    print_status "ERROR" "Neither wget nor curl found. Please install one of them."
    exit 1
fi

# Download each pyrho file
echo -e "\nDownloading pyrho LD block files..."
DOWNLOAD_SUCCESS=0
DOWNLOAD_TOTAL=${#PYRHO_FILES[@]}

for file in "${PYRHO_FILES[@]}"; do
    local_path="$DATA_DIR/$file"
    remote_url="$REPO_URL/$file"
    
    echo -e "\nDownloading $file..."
    
    # Check if file already exists
    if [ -f "$local_path" ]; then
        print_status "WARN" "File already exists: $local_path"
        read -p "Overwrite? (y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            print_status "OK" "Skipping $file"
            ((DOWNLOAD_SUCCESS++))
            continue
        fi
    fi
    
    # Download the file
    if [ "$DOWNLOAD_CMD" = "wget" ]; then
        if wget -q --show-progress -O "$local_path" "$remote_url"; then
            print_status "OK" "Downloaded: $file"
            ((DOWNLOAD_SUCCESS++))
        else
            print_status "ERROR" "Failed to download: $file"
            rm -f "$local_path"  # Remove partial download
        fi
    elif [ "$DOWNLOAD_CMD" = "curl" ]; then
        if curl -L --progress-bar -o "$local_path" "$remote_url"; then
            print_status "OK" "Downloaded: $file"
            ((DOWNLOAD_SUCCESS++))
        else
            print_status "ERROR" "Failed to download: $file"
            rm -f "$local_path"  # Remove partial download
        fi
    fi
done

# Verify downloads
echo -e "\n=============================================="
echo "Verifying downloaded files..."
echo "=============================================="

for file in "${PYRHO_FILES[@]}"; do
    local_path="$DATA_DIR/$file"
    if [ -f "$local_path" ]; then
        file_size=$(stat -c%s "$local_path" 2>/dev/null || echo "0")
        if [ "$file_size" -gt 0 ]; then
            print_status "OK" "$file (${file_size} bytes)"
        else
            print_status "ERROR" "$file is empty"
        fi
    else
        print_status "ERROR" "$file not found"
    fi
done

# Summary
echo -e "\n=============================================="
echo "Download Summary"
echo "=============================================="
echo "Downloaded: $DOWNLOAD_SUCCESS/$DOWNLOAD_TOTAL files"

if [ $DOWNLOAD_SUCCESS -eq $DOWNLOAD_TOTAL ]; then
    print_status "OK" "All pyrho LD block files downloaded successfully!"
    
    # Show usage information
    echo -e "\n${GREEN}Usage Information:${NC}"
    echo "The downloaded files are located in: $DATA_DIR/"
    echo ""
    echo "Available LD block files:"
    for file in "${PYRHO_FILES[@]}"; do
        echo "  - $file ($(echo $file | cut -d'_' -f2) population)"
    done
    echo ""
    echo "Update your config.yaml to point to these files:"
    echo "  ld_blocks: \"$DATA_DIR/pyrho_EUR_LD_blocks.bed\"  # for European ancestry"
    
    exit 0
else
    print_status "ERROR" "Some downloads failed. Please check your internet connection and try again."
    exit 1
fi