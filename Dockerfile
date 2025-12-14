FROM python:3.9-slim

# Set working directory inside container
WORKDIR /app

# 1. Install System Deps
RUN apt-get update && apt-get install -y \
    build-essential curl pkg-config libssl-dev git \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# 3. Install Python Deps
RUN pip install --no-cache-dir \
    pandas numpy scikit-learn scipy matplotlib seaborn joblib

# 4. Copy Source Code
# We copy files from host (left) to container (right)
COPY src/ ./src/
COPY Cargo.toml .
# Copy python scripts into a dedicated folder
COPY python/ ./python_scripts/

# 5. Build Rust
# We create a dummy build to handle the folder structure
RUN mkdir -p target/release && \
    cargo build --release

# 6. Setup Environment
# Move binary to global path
RUN mv target/release/gotchao_core /usr/local/bin/gotchao_core && \
    chmod +x /usr/local/bin/gotchao_core && \
    rm -rf target

# Set Python path so scripts find each other
ENV PYTHONPATH="/app/python_scripts:${PYTHONPATH}"

# 7. Entrypoint
# This makes the container executable like a script
ENTRYPOINT ["python3", "/app/python_scripts/run_GoTChAo.py"]