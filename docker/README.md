# Docker Usage Guide

This directory contains Docker configuration files for the FEM1D project.

## Quick Start

### Build and Run with Docker Compose

```bash
# Enter the docker/ directory
cd docker

# Build the images
docker-compose build --parallel

# Run TomlMain with default config
docker-compose run --rm fem-runtime

# Or pass a custom command, it will be parsed by entrypoint.sh
docker-compose run --rm fem-runtime help

# Run all good examples
docker-compose run --rm fem-test

# Interactive development environment
docker-compose run --rm -it fem-dev

# Sequential version for performance comparison
docker-compose run --rm fem-sequential

# Speedup analysis
docker-compose run --rm fem-speedup

# Interactive bash shell
docker-compose run --rm -it fem-bash
```

Note: 
- All runtime images will mount local `build/output` directory. 
- fem-runtime and fem-bash will mount local `config` directory as well.
- fem-dev will mount the whole project folder instead.

### Build and Run with Docker

```bash
# Build the images
docker build -t fem:runtime --target runtime .
docker build -t fem:dev --target development .

# Run with entrypoint commands
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime examples
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime speedup
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime help # List all commands

# Interactive development
docker run -it --rm -v $(pwd):/app fem:dev
```

## Available Commands (via entrypoint.sh)

The runtime image uses an entrypoint script that accepts these commands:

- **`toml [config_file]`** - Run TomlMain (default: config/good-examples/polynomial_1d.toml)
- **`sequential [config_file]`** - Run sequentialTomlMain (default: config/good-examples/polynomial_3d.toml)
- **`speedup [repeat] [config_dir] [output_file]`** - Run speedup analysis (default: repeat=3)
- **`examples`** - Run all good examples
- **`bash`** - Start interactive bash shell (use with `-it` flag)
- **`help`** - Show available commands

## Available Images

- **runtime**: Minimal image (~300MB) with runtime dependencies, built executables, and entrypoint script
- **development**: Full development environment (~1.3GB) with build tools, debuggers, and source code

## Docker Compose Services

- **fem-runtime**: Run TomlMain with polynomial_1d.toml config
- **fem-bash**: Interactive bash shell with runtime environment
- **fem-dev**: Interactive development environment with full build tools
- **fem-test**: Run all good examples automatically
- **fem-sequential**: Run sequential version with polynomial_3d.toml
- **fem-speedup**: Run speedup analysis with default parameters

## Volume Mounts

- `../output:/app/build/output` - Results and output files (runtime containers)
- `../config:/app/config` - Configuration files
- `..:/app` - Full source code (development only)
- `/app/build` - Anonymous volume for build artifacts (development)

## Usage Examples

### Running Different Commands

```bash
# Run specific TOML configurations
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime toml config/good-examples/trigonometric_2d.toml
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime sequential config/good-examples/gaussian_3d.toml

# Run speedup analysis with custom parameters
docker run --rm -v $(pwd)/output:/app/build/output fem:runtime speedup 5 ../config/speedup-analysis ../my_results.csv

# Interactive shell access
docker run --rm -it -v $(pwd)/output:/app/build/output fem:runtime bash

# Show help
docker run --rm fem:runtime help
```

### With Docker Compose

```bash
# Override default commands
docker-compose run --rm fem-runtime toml config/good-examples/exponential_3d.toml
docker-compose run --rm fem-runtime speedup 10
docker-compose run --rm fem-runtime help

# Run services with default commands
docker-compose run --rm fem-test     # Run all examples
docker-compose run --rm fem-speedup  # Speedup analysis
```

## Environment Variables

The entrypoint script uses these environment variables with defaults:

- **Default config files**: polynomial_1d.toml (toml), polynomial_3d.toml (sequential)
- **Speedup analysis**: REPEAT=3, config/speedup-analysis, speedup_results.csv
- **OpenMP**: Automatically uses all available CPU cores via `$(nproc)`

## Building for Different Architectures

```bash
# Build for ARM64 (Apple Silicon)
docker buildx build --platform linux/arm64 -t fem:arm64 --target runtime .

# Build for AMD64 (Intel/AMD)
docker buildx build --platform linux/amd64 -t fem:amd64 --target runtime .

# Multi-platform build
docker buildx build --platform linux/amd64,linux/arm64 -t fem:multi --target runtime .
```

## Development Workflow

1. Start development container:
   ```bash
   docker-compose run --rm -it fem-dev
   ```

2. Inside the container, rebuild when needed:
   ```bash
   cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)
   ```

3. Run tests:
   ```bash
   ./17_fem1d_17_fem1d_test
   ```

4. Run specific configurations:
   ```bash
   ./TomlMain config/good-examples/polynomial_3d.toml
   ./sequentialTomlMain config/good-examples/trigonometric_2d.toml
   ```

## Performance Notes

- **Runtime image**: Optimized for size (~72MB vs ~1.3GB for development)
- **OpenMP enabled**: Parallel execution with `libgomp1`
- **Multi-core builds**: Uses `make -j$(nproc)` for faster compilation
- **Sequential comparison**: Use `sequentialTomlMain` for single-threaded benchmarks
- **Output organization**: Results saved to `output/good_examples/` directory

## File Structure

```
/app/
├── build/
│   ├── TomlMain              # Parallel FEM solver
│   ├── sequentialTomlMain    # Sequential FEM solver
│   └── output/               # Results output directory
├── config/                   # TOML configuration files
├── mesh/                     # Mesh files (.msh)
├── speedup_analysis/         # Performance analysis scripts
├── run_good_examples.sh      # Script to run all examples
└── entrypoint.sh            # Docker entrypoint with commands
```

## Troubleshooting

### Common Issues

**Build failures:**
- Check that all source files are present and not excluded in `.dockerignore`
- Ensure CMake version >=3.26 (automatically installed via pip in container)
- Verify all dependencies are available in Ubuntu 22.04 repositories

**Permission issues:**
- Runtime containers run as `femuser` (non-root)
- Development containers run as `devuser` (non-root)
- Ensure host output directory has write permissions: `chmod 755 ./output`

**TTY warnings:**
- Use `-it` flags for interactive containers: `docker run -it fem:runtime bash`
- TTY environment variables are pre-configured in the images

**Volume mount issues:**
- Runtime working directory is `/app/build/`
- Output should be mounted to `/app/build/output`
- Config files should be mounted to `/app/config`

**Memory/performance:**
- Container uses all available CPU cores by default
- Limit with `--cpus="4"` if needed
- For large problems, ensure sufficient Docker memory allocation

### Getting Help

```bash
# Show available entrypoint commands
docker run --rm fem:runtime help

# Check container status
docker-compose ps

# View container logs
docker-compose logs fem-runtime

# Enter running development container
docker-compose exec fem-dev bash
```

### Mesh File Errors

If you encounter mesh file errors:
- Verify mesh files are in `/app/mesh/` inside container
- Check mesh file format compatibility (Gmsh v2.2 format)
- Ensure mesh files are not excluded in `.dockerignore`

### Output Files

All output files are organized as follows:
- **Good examples**: `output/good_examples/[problem_name].*`
- **Speedup analysis**: `output/speedup_results.csv`
- **Custom runs**: `output/[custom_name].*`

The output directory is automatically created and mounted from the host system.
