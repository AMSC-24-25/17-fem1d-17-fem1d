# Multi-stage Dockerfile for 17-fem1d-17-fem1d project
FROM ubuntu:22.04 AS builder

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    pkg-config \
    libeigen3-dev \
    libgomp1 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install newer CMake from pip (Ubuntu 22.04 has older CMake)
RUN pip3 install cmake>=3.26

# Set working directory
WORKDIR /app

# Copy source code
COPY .. .

# Create build directory and configure
RUN mkdir -p build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release

# Build the project
RUN cd build && make -j$(nproc) TomlMain && make -j$(nproc) sequentialTomlMain

# Runtime stage - smaller image for running
FROM ubuntu:22.04 AS runtime

# Set terminal environment variables to prevent TTY warnings
ENV DEBIAN_FRONTEND=
ENV TERM=xterm-256color
ENV COLUMNS=80
ENV LINES=24

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    libeigen3-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN useradd -m -s /bin/bash femuser

# Set working directory
WORKDIR /app/
RUN mkdir -p build

# Copy built executables and necessary files from builder stage
COPY --from=builder /app/build/TomlMain /app/build/
COPY --from=builder /app/build/sequentialTomlMain /app/build/
COPY --from=builder /app/config /app/config/
COPY --from=builder /app/mesh /app/mesh/
COPY --from=builder /app/run_good_examples.sh /app/
# Copy necessary files for running speedup analysis
COPY --from=builder /app/speedup_analysis /app/speedup_analysis/
# Copy entrypoint script
COPY entrypoint.sh /app/

# Create output directory
RUN mkdir -p build/output && \
    chown -R femuser:femuser /app

# Make scripts executable
RUN chmod +x ./*.sh
RUN chmod +x ./speedup_analysis/*.sh

# Switch to non-root user
USER femuser

# Set default command and work directory
WORKDIR /app/build
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["speedup"]

# Development stage - includes all build tools
FROM builder AS development

# Set terminal environment variables to prevent TTY warnings  
ENV TERM=xterm-256color
ENV COLUMNS=80
ENV LINES=24

# Install additional development tools
RUN apt-get update && apt-get install -y \
    gdb \
    valgrind \
    vim \
    nano \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for development
RUN useradd -m -s /bin/bash devuser && \
    chown -R devuser:devuser /app

USER devuser

WORKDIR /app

# Default command for development
CMD ["/bin/bash"]
