#!/bin/bash

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Run this script from scripts directory.
if [[ "$(basename $(pwd))" != "scripts" ]]; then
    echo -e "${RED}This script should be run from the scripts directory.${NC}"
    exit 1
fi

cd ..
rm -rf build
mkdir build
cd build
cmake ..
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Build successful!${NC}"
else
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi