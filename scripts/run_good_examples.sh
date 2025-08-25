#!/bin/bash

# Script to run all TOML configuration files in config/good-examples folder
# This script executes the FEM solver on all test problems with known exact solutions

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
EXECUTABLE="./build/TomlMain"
CONFIG_DIR="./config/good-examples"
OUTPUT_DIR="./build/output/good_examples"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable $EXECUTABLE not found!${NC}"
    echo "Please build the project first with 'make' in the build directory."
    exit 1
fi

# Check if config directory exists
if [ ! -d "$CONFIG_DIR" ]; then
    echo -e "${RED}Error: Config directory $CONFIG_DIR not found!${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Running FEM solver on good examples${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Initialize counters
total_files=0
successful_runs=0
failed_runs=0

# Process all .toml files in the good-examples directory
for config_file in "$CONFIG_DIR"/*.toml; do
    if [ -f "$config_file" ]; then
        total_files=$((total_files + 1))
        filename=$(basename "$config_file")
        name_without_ext="${filename%.toml}"
        
        echo -e "${YELLOW}Running: $filename${NC}"
        echo "Config: $config_file"
        
        # Run the FEM solver
        if "$EXECUTABLE" "$config_file" > "$OUTPUT_DIR/${name_without_ext}.log" 2>&1; then
            echo -e "${GREEN}✓ Success${NC}"
            successful_runs=$((successful_runs + 1))
        else
            echo -e "${RED}✗ Failed${NC}"
            failed_runs=$((failed_runs + 1))
            echo "  Check log: $OUTPUT_DIR/${name_without_ext}.log"
        fi
        
        echo ""
    fi
done

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo "Total files processed: $total_files"
echo -e "Successful runs: ${GREEN}$successful_runs${NC}"
echo -e "Failed runs: ${RED}$failed_runs${NC}"
echo ""

mkdir -p $OUTPUT_DIR/csv
mkdir -p $OUTPUT_DIR/log
mv $OUTPUT_DIR/*.log $OUTPUT_DIR/log/
mv $OUTPUT_DIR/*.csv $OUTPUT_DIR/csv/

echo "Output vtu saved to: $OUTPUT_DIR"
echo "Output csv saved to: $OUTPUT_DIR/csv"
echo "Output log saved to: $OUTPUT_DIR/log"
echo ""

if [ $failed_runs -eq 0 ]; then
    echo -e "${GREEN}All tests completed successfully!${NC}"
    exit 0
else
    echo -e "${YELLOW}Some tests failed. Check the log files for details.${NC}"
    exit 1
fi
