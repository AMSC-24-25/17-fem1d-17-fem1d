#!/bin/bash

# Docker entrypoint script for multiple command options
# If running manually (e.g. not in Docker), run from build directory
set -e

case "$1" in
    "toml")
        shift
        CONFIG_FILE="${1:-config/good-examples/polynomial_1d.toml}"
        echo "Running TomlMain with config: $CONFIG_FILE"
        exec ./TomlMain "$CONFIG_FILE"
        ;;
    "sequential")
        shift
        CONFIG_FILE="${1:-config/good-examples/polynomial_3d.toml}"
        echo "Running sequentialTomlMain with config: $CONFIG_FILE"
        exec ./sequentialTomlMain "$CONFIG_FILE"
        ;;
    "speedup")
        REPEAT="${2:-3}"
        CONFIG_DIR="${3:-../config/speedup-analysis}"
        OUTPUT_FILE="${4:-../speedup_results.csv}"
        echo "Running speedup analysis (REPEAT=$REPEAT)"
        exec bash ../speedup_analysis/run_speedup.sh ./TomlMain "$CONFIG_DIR" "$OUTPUT_FILE" ./sequentialTomlMain
        ;;
    "examples")
        echo "Running all good examples"
        exec ../run_good_examples.sh
        ;;
    "bash")
        echo "Starting interactive bash shell"
        echo "If this closes immediately, try running with docker flag '-it'"
        exec /bin/bash
        ;;
    "help"|"--help"|"-h")
        echo "Available commands:"
        echo "  toml [config_file]           - Run TomlMain (default: config/good-examples/polynomial_1d.toml)"
        echo "  sequential [config_file]     - Run sequentialTomlMain (default: config/good-examples/polynomial_3d.toml)"
        echo "  speedup [repeat] [config_dir] [output_file] - Run speedup analysis (default: repeat=3)"
        echo "  examples                     - Run all good examples"
        echo "  bash                         - Start interactive bash shell (remember -it flag)"
        echo "  help                         - Show this help"
        ;;
    *)
        echo "Unknown command: $1"
        echo "Use 'help' to see available commands"
        exit 1
        ;;
esac
