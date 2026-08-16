# Dockerfile for PSI (Parallel Spatial Indexes) Artifact Evaluation
FROM ubuntu:24.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Build essentials
    build-essential \
    cmake \
    g++-14 \
    gcc-14 \
    git \
    wget \
    # Boost library
    libboost-all-dev \
    # Optional performance libraries
    libjemalloc-dev \
    numactl \
    # Python dependencies
    python3 \
    python3-pip \
    python3-venv \
    # R for plotting
    r-base \
    r-base-dev \
    # Utilities
    time \
    vim \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Install R packages for plotting
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'ggsci', 'ggpubr', 'latex2exp', 'scales'), repos='http://cran.rstudio.com/')"

# Install Python packages
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    pandas \
    matplotlib

# Set C++ compiler version
ENV CXX=g++-14
ENV CC=gcc-14

# Copy the project files
COPY . /workspace/SpaceTreeLib

# Set working directory to project
WORKDIR /workspace/SpaceTreeLib

# Initialize git submodules
RUN git submodule update --init --recursive || true

# Binaries baked into a distributed image must run on machines other than the
# one that built them, so the image is compiled portably by default. That also
# means these binaries are NOT tuned for the host: recompile inside the
# container with -DPSI_NATIVE_ARCH=ON before any timing run, or build the
# image locally with --build-arg PSI_NATIVE_ARCH=ON.
ARG PSI_NATIVE_ARCH=OFF

# Create build directory and compile the project
RUN mkdir -p build && \
    cd build && \
    cmake -DDEBUG=OFF -DCGAL=OFF -DJEMA=OFF \
          -DPSI_NATIVE_ARCH=${PSI_NATIVE_ARCH} .. && \
    make -j$(nproc)

# Create directories for data and results
RUN mkdir -p script_ae/logs && \
    mkdir -p script_ae/data && \
    mkdir -p script_ae/plots

# Set default environment variables
ENV DATA_PREFIX=/data
ENV NODE_SIZE=1000000000

# Set the entry point to bash for interactive use
CMD ["/bin/bash"]
