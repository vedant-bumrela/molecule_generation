# Quick Start Guide

This molecular docking system can be run in three ways:

## Option 1: Docker (Recommended - Easiest)

The easiest way to run the complete system with web interface and batch processing:

```bash
# Build and start all services (Flask API, Redis, Worker)
docker-compose up -d

# Access the web interface
# Open browser to: http://localhost:5000

# View logs
docker-compose logs -f app

# Stop services
docker-compose down
```

**Requirements:**
- Docker installed
- Docker Compose installed
- NVIDIA Docker runtime (for GPU support with ML models)
- NVIDIA GPU with CUDA support

## Option 2: Web API (Local Development)

Run the Flask web API with web interface:

```bash
# Install dependencies (if not already done)
pip install -r requirements.txt

# You'll also need system dependencies:
# - AutoDock Vina
# - Open Babel
# - RDKit (from conda-forge)

# Run the Flask API
python -m src.api.app

# Or run with gunicorn for production
gunicorn -w 4 -b 0.0.0.0:5000 src.api.app:app

# Access: http://localhost:5000
```

**In a separate terminal**, run the worker to process jobs:

```bash
# Start the worker
python -m src.api.worker

# Make sure Redis is running locally
# On Ubuntu/Debian: sudo apt-get install redis-server && redis-server
# On MacOS: brew install redis && redis-server
# On Windows: Download from https://redis.io/download
```

## Option 3: Command Line Tool (CLI)

Run single docking jobs from the command line:

```bash
# Basic usage with protein and ligand files
python dock.py --protein path/to/protein.pdb --ligand path/to/ligand.sdf --out ./results

# Or use SMILES string
python dock.py --protein path/to/protein.pdb --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --out ./results

# With specific binding site
python dock.py \
  --protein path/to/protein.pdb \
  --ligand path/to/ligand.sdf \
  --center_x 10.0 --center_y 20.0 --center_z 15.0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --exhaustiveness 16 \
  --num-modes 20 \
  --out ./results

# See all options
python dock.py --help
```

## Prerequisites

### System Dependencies

**All platforms need:**
- Python 3.8+
- AutoDock Vina (download from https://vina.scripps.edu/)
- Open Babel (https://openbabel.org/docs/dev/Installation/install.html)

**For Windows:**
```powershell
# Download AutoDock Vina from https://vina.scripps.edu/download/vina
# Extract and add to PATH

# Install Open Babel
choco install openbabel
```

**For Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git wget curl
sudo apt-get install -y libboost-all-dev libopenbabel-dev python3-dev python3-pip

# Install AutoDock Vina
wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
tar -xzf autodock_vina_1_1_2_linux_x86.tgz
sudo mv autodock_vina_1_1_2_linux_x86/bin/* /usr/local/bin/
```

**For MacOS:**
```bash
# Install via Homebrew
brew install open-babel
brew install autodock-vina
```

### Python Dependencies

**Using pip:**
```bash
pip install -r requirements.txt
```

**Using conda (recommended for RDKit):**
```bash
conda env create -f environment.yml
conda activate docking
```

### ML Model Dependencies (Optional)

The ML models (GNINA, EquiBind, DiffDock) will be automatically cloned in Docker but need manual setup for local runs:

```bash
# GNINA will be auto-cloned in Docker
# For local, you can skip or install manually

# EquiBind and DiffDock
git clone https://github.com/HannesStark/EquiBind.git
git clone https://github.com/gcorso/DiffDock.git

# Download their pretrained models separately
```

## Usage Examples

### Web Interface

1. Start the system: `docker-compose up -d`
2. Open http://localhost:5000
3. Upload protein PDB file
4. Upload ligand SDF or enter SMILES
5. Select docking method
6. Submit job
7. View results in 3D viewer

### Batch Processing

1. Go to "Batch Processing" tab in web interface
2. Upload protein structure
3. Upload CSV file with SMILES column or SDF library
4. Select parameters
5. Submit batch job
6. Monitor progress and download results

### Command Line

```bash
# Simple docking
python dock.py --protein 1abc.pdb --ligand ligand.sdf

# With SMILES
python dock.py --protein 1abc.pdb --smiles "CCO" 

# High-throughput with custom parameters
python dock.py \
  --protein protein.pdb \
  --smiles "CCO" \
  --center_x 0.0 --center_y 0.0 --center_z 0.0 \
  --exhaustiveness 32 \
  --num-modes 20 \
  --out ./docking_results
```

## Troubleshooting

### Redis Connection Error
```bash
# Make sure Redis is running
redis-cli ping  # Should return PONG

# Or in Docker, check logs
docker-compose logs redis
```

### Vina Not Found
```bash
# Check if vina is in PATH
which vina  # Linux/Mac
where vina  # Windows

# If not found, add to PATH or specify in config
```

### CUDA/GPU Issues
```bash
# Check NVIDIA Docker runtime
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi

# If this fails, install NVIDIA Docker runtime
```

### Import Errors
```bash
# Make sure you're in the correct directory
cd c:\docking

# Activate conda environment if using
conda activate docking

# Reinstall dependencies
pip install -r requirements.txt
```

## Getting Help

- Check logs: `docker-compose logs -f`
- Enable verbose mode: `--verbose` flag
- API health check: `curl http://localhost:5000/api/health`

