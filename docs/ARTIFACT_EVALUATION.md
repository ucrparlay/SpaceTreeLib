## Docker Setup (Recommended for Artifact Evaluation)

**Hardware Requirements:**
- **Memory**: at least 512GB RAM (preferable >1T) for full experiments
- **CPU**: All available cores (preferable >128 cores)
- **Disk**: 100GB+ free space

**Warning**: The scripts for results merging and figures plotting may not execute properly if a baseline failed to generate the log and exit correctly, i.e., the resource requested exceeds the RAM. Please make sure the machine where the artifact is evaluated satisfied above requirements. Thanks. 

**Before timing anything:** the published image is compiled portably, so its
binaries run on any x86-64 machine but are *not* tuned for yours. Recompile
inside the container before a measurement run:

```bash
cmake -S . -B build -DDEBUG=OFF -DCGAL=OFF -DPSI_NATIVE_ARCH=ON && \
  cmake --build build -j
```

Building the image locally with `--build-arg PSI_NATIVE_ARCH=ON` does the same
thing. This matters: the numbers in the papers were measured with
`-march=native -mcx16`, and a portable build will be slower.

**Quick Start:**

```bash
# Pull the prebuilt image
./docker-run.sh pull

# Run customized experiments (specify where to store generated data and the synthetic benchmark size)
./docker-run.sh run --data-path /mnt/large-disk/data --node-size 1000000

# Or run the full evaluation in paper
./docker-run.sh full --data-path /mnt/large-disk/data
```

**Parameters:**

- `--data-path PATH` - Host directory for generated datasets (required)
- `--download-real-world-data` - Download real-world data on host (default: true)
- `--no-download-real-world-data` - Skip downloading real-world data
- `--node-size SIZE` - Number of data points (default: 1'000'000'000, aka. 1 billion)
- `--memory SIZE` - Memory limit (e.g., 512g), generally you don't want to explicitly set this value
- `--cpus NUM` - CPU limit (default: all cores), generally you don't want set it as well.

**Results:** Found in `script_ae/data/`, `script_ae/logs/`, `script_ae/plots/` directories on the host.

**Documentation:** See [docker reference](https://github.com/ucrparlay/SpaceTreeLib/blob/main/docs/DOCKER_QUICK_REFERENCE.md) for complete guides to use docker. More introduction about the project can be found in the [Github page](https://github.com/ucrparlay/SpaceTreeLib).
