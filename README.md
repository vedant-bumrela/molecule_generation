# AI-Enabled Molecular Docking System

A comprehensive platform for molecular docking with AI-powered pose prediction and scoring.

## Features

- **Multiple Docking Methods**:
  - Classical docking with AutoDock Vina
  - ML-based rescoring with GNINA
  - ML-based pose prediction with EquiBind and DiffDock

- **Batch Processing**:
  - High-throughput screening (HTS) capabilities
  - Process large compound libraries (SMILES or SDF)
  - Distributed task processing with Redis queue

- **Web Interface**:
  - Interactive 3D visualization with NGL Viewer
  - Job submission and management
  - Results visualization and download

## Architecture

- **Backend**: Flask API with Redis task queue
- **Frontend**: HTML/CSS/JS with Bootstrap and NGL Viewer
- **Processing**: Multi-threaded batch processing
- **Containerization**: Docker with CUDA support

## Installation

### Using Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/molecular-docking.git
cd molecular-docking

# Build and start the containers
docker-compose up -d
```

### Manual Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/molecular-docking.git
cd molecular-docking

# Create and activate a conda environment
conda env create -f environment.yml
conda activate docking

# Install dependencies
pip install -r requirements.txt

# Start the application
python src/api/app.py
```

## Usage

### Web Interface

Access the web interface at http://localhost:5000

### Single Docking

1. Upload a protein structure (PDB format)
2. Upload a ligand (SDF format) or enter SMILES
3. Select docking method
4. Set docking parameters (optional)
5. Submit job and view results

### Batch Processing

1. Upload a protein structure (PDB format)
2. Upload a compound library (CSV with SMILES or SDF)
3. Select docking method
4. Set number of worker threads
5. Submit batch job and monitor progress

## API Endpoints

- `GET /`: Web interface
- `POST /submit`: Submit single docking job
- `POST /submit/batch`: Submit batch docking job
- `GET /status/<job_id>`: Get job status
- `GET /batch/results/<job_id>`: Get batch job results
- `GET /download/<job_id>/<filename>`: Download result files
- `GET /jobs`: List all jobs

## Dependencies

- AutoDock Vina
- Open Babel
- RDKit
- PyTorch
- Flask
- Redis
- GNINA (optional)
- EquiBind (optional)
- DiffDock (optional)

## License

MIT

## Citation

If you use this software in your research, please cite:

```
@software{molecular_docking,
  author = {Your Name},
  title = {AI-Enabled Molecular Docking System},
  year = {2023},
  url = {https://github.com/yourusername/molecular-docking}
}
```