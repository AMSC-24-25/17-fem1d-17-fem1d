#!/bin/bash

# Run this script from root, after compiling.
echo "This script should be run from root, after compiling."

echo "Testing TOML configuration files..."

# Create output directory if it doesn't exist
OUTPUT_DIR="build/output"
mkdir -p $OUTPUT_DIR

# Array of config files to test
configs=("debug_test" "test_simple" "poisson" "problem_1d" "problem_2d" "transport_dominated")

echo "Testing configurations:"
for config in "${configs[@]}"; do
    echo "----------------------------------------"
    echo "Testing: $config.toml"
    echo "----------------------------------------"
    
    ./build/TomlMain ./config/$config.toml
    
    if [ $? -eq 0 ]; then
        echo "✓ $config.toml - SUCCESS"
        # Check if output files were created
        if [ -f "$OUTPUT_DIR/${config}.csv" ]; then
            echo "  - CSV file created: $OUTPUT_DIR/${config}.csv"
        fi
        if [ -f "$OUTPUT_DIR/${config}.vtu" ]; then
            echo "  - VTU file created: $OUTPUT_DIR/${config}.vtu"
        fi
    else
        echo "✗ $config.toml - FAILED"
    fi
    echo ""
done

echo "Testing completed!"
echo "Check the $OUTPUT_DIR directory for generated files."
