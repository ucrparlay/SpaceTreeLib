# Docker Quick Reference

## Build and Run

```bash
# 1. Build image
./docker-run.sh build

# 2. Run experiments (requires --data-path)
./docker-run.sh run --data-path /mnt/large-disk/data

# 3. Full evaluation
./docker-run.sh full --data-path /mnt/large-disk/data
```

## Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `--data-path PATH` | Host directory for datasets (required) | - |
| `--node-size SIZE` | Number of data points | 1000000000 |
| `--memory SIZE` | Memory limit (e.g., 512g) | No limit |
| `--cpus NUM` | CPU cores limit | All cores |

## Examples

```bash
# Full experiments with 512GB RAM
./docker-run.sh full --data-path /mnt/ssd/data --memory 512g

# Custom dataset size
./docker-run.sh run --data-path /data --node-size 100000

# Interactive mode
./docker-run.sh run --data-path /scratch/exp
# Then inside: ./run.sh /data 1000000000
```

## Results Location

| Type | Location |
|------|----------|
| Generated data | Host path specified in `--data-path` |
| Summarized results | `./data/` |
| Logs | `./logs/` |
| Plots | `./plots/` |

## Using Docker Compose

```bash
# Edit .env file
cp .env.example .env
# Set DATA_MOUNT_PATH=/your/data/path

# Start
docker-compose up -d

# Enter container
docker-compose exec psi-ae bash

# Stop
docker-compose down
```

## Cleanup

```bash
# Remove image
./docker-run.sh clean

# Remove generated data
rm -rf /your/data/path/*
rm -rf results/ logs/ plots/
```

## Hardware Requirements

| Experiment Size | RAM | CPU | Disk |
|----------------|-----|-----|------|
| Full (1B points) | 512GB+ | All cores | 100GB+ |

## Troubleshooting

**Issue:** `numactl: Operation not permitted`  
**Fixed:** Container runs with `--cap-add=SYS_NICE`

**Issue:** Files not visible on host  
**Solution:** Check mounted directories (data/, logs/, plots/)

**Issue:** Out of memory  
**Solution:** Reduce `--node-size` or increase system RAM

**Issue:** Script changes not reflected  
**Solution:** Rebuild image with `./docker-run.sh build`
