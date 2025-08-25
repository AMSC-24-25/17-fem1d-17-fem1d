#!/bin/bash

# Run this script from root, after compiling.
if [[ "$(basename $(pwd))" != "scripts" ]]; then
    echo -e "${RED}This script should be run from the scripts directory, after compiling.${NC}"
    exit 1
fi
cd ../build

echo "Testing TOML configuration files..."

# Create output directory if it doesn't exist
OUTPUT_DIR="output"
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/log
mkdir -p $OUTPUT_DIR/csv

# Array of config files to test - automatically discover all .toml files
configs=()
while IFS= read -r -d '' file; do
    filename=$(basename "$file" .toml)
    configs+=("$filename")
done < <(find ../config -name "*.toml" -type f -print0)

echo "Found ${#configs[@]} configuration files in ./config/"
echo ${configs[@]}

# Initialize counters
total_tests=0
success_count=0
failure_count=0

echo "Testing configurations:"
for config in "${configs[@]}"; do
    config_file="../config/$config.toml"
    echo "----------------------------------------"
    echo "Testing: $config_file"
    echo "----------------------------------------"

    total_tests=$((total_tests + 1))
    ./TomlMain "$config_file" > "./$OUTPUT_DIR/log/${config}.log" 2>&1

    if [ $? -eq 0 ]; then
        echo "✓ $config.toml - SUCCESS"
        success_count=$((success_count + 1))
        # Check if output files were created
        if [ -f "$OUTPUT_DIR/${config}.csv" ]; then
            mv "$OUTPUT_DIR/${config}.csv" "$OUTPUT_DIR/csv/${config}.csv"
            echo "  - CSV file created: $OUTPUT_DIR/csv/${config}.csv"
        fi
        if [ -f "$OUTPUT_DIR/${config}.vtu" ]; then
            echo "  - VTU file created: $OUTPUT_DIR/${config}.vtu"
        fi
    else
        echo "✗ $config.toml - FAILED"
        failure_count=$((failure_count + 1))
        echo "  - Check log: build/$OUTPUT_DIR/log/${config}.log"
        echo "  - Last log line: $(tail -n 1 ./$OUTPUT_DIR/log/${config}.log)"
    fi
    echo ""
done

echo "========================================="
echo "Testing Summary"
echo "========================================="
echo "Total tests run: $total_tests"
echo "Successful: $success_count"
echo "Failed: $failure_count"
echo ""

if [ $failure_count -eq 0 ]; then
    echo "All tests passed!"
    exit_code=0
else
    echo "Some tests failed. Check the log files in $OUTPUT_DIR/log/ for details."
    exit_code=1
fi

echo "\nTesting completed!"
echo "Check the $OUTPUT_DIR directory for generated files."

exit $exit_code
