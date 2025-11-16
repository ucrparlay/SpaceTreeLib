## Docker Setup (Recommended for Artifact Evaluation)

**Hardware Requirements:**
- **Memory**: at least 512GB RAM (preferable >1T) for full experiments
- **CPU**: All available cores (preferable >128 cores)
- **Disk**: 100GB+ free space

**Quick Start:**
```bash
# Build the Docker image
./docker-run.sh build

# Run experiments (specify where to store generated data)
./docker-run.sh run --data-path /mnt/large-disk/data

# Or run full evaluation
./docker-run.sh full --data-path /mnt/large-disk/data
```

**Parameters:**
- `--data-path PATH` - Host directory for generated datasets (required)
- `--node-size SIZE` - Number of data points (default: 1'000'000'000, aka. 1 billion)
- `--memory SIZE` - Memory limit (e.g., 512g), generally you don't want to explicitly set this value
- `--cpus NUM` - CPU limit (default: all cores), generally you don't want set it as well.

**Results:** Found in `data/`, `logs/`, `plots/` directories on host.

**Documentation:** See [doc/](doc/) folder for complete guides.
