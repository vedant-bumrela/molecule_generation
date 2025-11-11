FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    libboost-all-dev \
    libopenbabel-dev \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install AutoDock Vina
RUN wget --no-check-certificate https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz \
    && tar -xzf autodock_vina_1_1_2_linux_x86.tgz \
    && mv autodock_vina_1_1_2_linux_x86/bin/vina /usr/local/bin/ \
    && mv autodock_vina_1_1_2_linux_x86/bin/vina_split /usr/local/bin/ \
    && rm -rf autodock_vina_1_1_2_linux_x86 autodock_vina_1_1_2_linux_x86.tgz

# Install Open Babel
RUN apt-get update && apt-get install -y --no-install-recommends \
    openbabel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Install PyTorch with CUDA support (with --no-deps to avoid hash verification)
RUN pip3 install --no-cache-dir --no-deps torch==1.13.1+cu117 torchvision==0.14.1+cu117 -f https://download.pytorch.org/whl/torch_stable.html

# Provide typing_extensions required by torch at import-time for some builds
RUN pip3 install --no-cache-dir typing_extensions==4.7.1

# Install PyTorch Geometric with compatible version
RUN pip3 install --no-cache-dir --no-deps torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.13.1+cu117.html

# Install dependencies that might have been skipped
RUN pip3 install --no-cache-dir numpy pandas scipy matplotlib scikit-learn

# Copy application code
COPY . .

# Clone ML repositories
RUN mkdir -p /app/ml_models && \
    cd /app/ml_models && \
    git clone https://github.com/gnina/gnina.git && \
    git clone https://github.com/HannesStark/EquiBind.git && \
    git clone https://github.com/gcorso/DiffDock.git

# Create directories for data and results
RUN mkdir -p /app/data/uploads /app/data/results

# Expose port
EXPOSE 5000

# Run the application as a module so relative imports work
CMD ["python3", "-m", "src.api.app"]