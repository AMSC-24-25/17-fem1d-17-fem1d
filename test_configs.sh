#!/bin/bash
echo "Testing TOML configuration files..."

# Create output directory if it doesn't exist
mkdir -p output

# Array of config files to test
configs=("debug_test" "test_simple" "poisson" "problem_1d" "problem_2d" "transport_dominated")

echo "Building project..."
cmake --build build
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo "Testing configurations:"
for config in "${configs[@]}"; do
    echo "----------------------------------------"
    echo "Testing: $config.toml"
    echo "----------------------------------------"
    
    ./build/TomlMain config/$config.toml
    
    if [ $? -eq 0 ]; then
        echo "✓ $config.toml - SUCCESS"
        # Check if output files were created
        if [ -f "output/${config}.csv" ]; then
            echo "  - CSV file created: output/${config}.csv"
        fi
        if [ -f "output/${config}.vtu" ]; then
            echo "  - VTU file created: output/${config}.vtu"
        fi
    else
        echo "✗ $config.toml - FAILED"
    fi
    echo ""
done

echo "Testing completed!"
echo "Check the output/ directory for generated files."
