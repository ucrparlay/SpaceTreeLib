# Docker Quick Reference

## Quick Start - Using Pre-built Image (Recommended)

The easiest way to use PSI is to pull the pre-built Docker image from GitHub Container Registry:

```bash
# 1. Pull the pre-built image (no build needed!)
./docker-run.sh pull

# 2. Run experiments (requires --data-path)
./docker-run.sh run --data-path /mnt/large-disk/data --node-size 1000000

# 3. Full evaluation
./docker-run.sh full --data-path /mnt/large-disk/data
```

## Alternative - Build Locally

If you want to modify the code or can't access the pre-built image:

```bash
# 1. Build image locally
./docker-run.sh build

# 2. Run experiments
./docker-run.sh run --data-path /mnt/large-disk/data
```

## Manual Docker Commands

You can also use Docker directly without the helper script:

```bash
# Pull the image
docker pull ghcr.io/ucrparlay/spacetreelib:latest

# Run interactively
docker run -it --rm \
  -v /ssd1:/ssd1 \
  -e DATA_PREFIX=/ssd1/ \
  -e NODE_SIZE=1000000 \
  ghcr.io/ucrparlay/spacetreelib:latest

# Run with custom settings
docker run -it --rm \
  --cap-add=SYS_NICE \
  -v $(pwd)/script_ae/logs:/workspace/SpaceTreeLib/script_ae/logs \
  -e DATA_PREFIX=/data \
  -e NODE_SIZE=1000000000 \
  ghcr.io/ucrparlay/spacetreelib:latest
```

## Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `--data-path PATH` | Host directory for datasets (required) | - |
| `--download-real-world-data` | Download real-world data on host | `true` |
| `--no-download-real-world-data` | Skip downloading real-world data | - |
| `--node-size SIZE` | Number of data points | 1000000000 |
| `--memory SIZE` | Memory limit (e.g., 512g) | No limit |
| `--cpus NUM` | CPU cores limit | All cores |

## Examples

```bash
# Full experiments
./docker-run.sh full --data-path /mnt/ssd/data

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
# Remove container, image, and generated files in script_ae/
./docker-run.sh clean

# Remove container, image, generated files AND real-world data in --data-path
./docker-run.sh clean --data-path /your/data/path

# Or manually remove the image
docker rmi ghcr.io/ucrparlay/spacetreelib:latest
```

## Hardware Requirements

| Experiment Size | RAM | CPU | Disk |
|----------------|-----|-----|------|
| Full (1B points) | 512GB+ | All cores | 100GB+ |

## Troubleshooting

**Issue:** `numactl: Operation not permitted`  
**Fixed:** Container runs with `--cap-add=SYS_NICE`

**Issue:** Files not visible on host  
**Solution:** Check mounted directories (script_ae/data/, script_ae/logs/, script_ae/plots/)

**Issue:** Out of memory  
**Solution:** Reduce `--node-size` or increase system RAM

**Issue:** Script changes not reflected  
**Solution:** Rebuild image with `./docker-run.sh build` or pull latest with `./docker-run.sh pull`

**Issue:** Cannot pull image (permission denied)  
**Solution:** The package may be private. Ask repository admin to make it public in GitHub package settings, or authenticate with:
```bash
echo YOUR_GITHUB_TOKEN | docker login ghcr.io -u YOUR_USERNAME --password-stdin
```
