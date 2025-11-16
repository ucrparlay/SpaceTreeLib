## Docker Setup (Recommended for Artifact Evaluation)

**Hardware Requirements:**
- **Memory**: at least 512GB RAM (preferable >1T) for full experiments
- **CPU**: All available cores (preferable >128 cores)
- **Disk**: 100GB+ free space

**Quick Start:**
```bash
# Build the Docker image
./docker-run.sh pull

# Run experiments (specify where to store generated data and the synthetic benchmark size)
./docker-run.sh run --data-path /mnt/large-disk/data --node-size 1000000

# Or run full evaluation
./docker-run.sh full --data-path /mnt/large-disk/data
```

**Parameters:**
- `--data-path PATH` - Host directory for generated datasets (required)
- `--node-size SIZE` - Number of data points (default: 1'000'000'000, aka. 1 billion)
- `--memory SIZE` - Memory limit (e.g., 512g), generally you don't want to explicitly set this value
- `--cpus NUM` - CPU limit (default: all cores), generally you don't want set it as well.

**Results:** Found in `script_ae/data/`, `script_ae/logs/`, `script_ae/plots/` directories on host.

**Documentation:** See [docker reference](docs/DOCKER_QUICK_REFERENCE.md) for complete guides to use docker.
