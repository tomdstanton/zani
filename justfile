# Zani Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-c"]

# Show available commands
default:
    @just --list

# Sync python dependencies and create the virtual environment using `uv`
sync:
    uv sync

# Build the Rust CLI and the Python extension module (development)
build: sync
    cargo build --workspace
    uv run maturin develop

# Build the Rust CLI and the Python extension module (release)
release: sync
    cargo build --workspace --release
    uv run maturin develop --release

# Generate CLI documentation from the binary help output
doc-cli: build
    mkdir -p docs
    echo "# CLI Reference" > docs/cli.md
    echo "\`\`\`text" >> docs/cli.md
    cargo run --bin zani -- --help >> docs/cli.md
    echo "\`\`\`" >> docs/cli.md

# Generate benchmark documentation
doc-bench: build
    mkdir -p docs
    echo "# Benchmarks" > docs/benchmarks.md
    echo "\`\`\`text" >> docs/benchmarks.md
    cargo bench --workspace >> docs/benchmarks.md || true
    echo "\`\`\`" >> docs/benchmarks.md

# Build all documentation and generate Zensical site
docs: build doc-cli doc-bench
    cargo doc --no-deps --workspace
    mkdir -p docs/rust-api
    cp -r target/doc/* docs/rust-api/
    cp README.md docs/index.md
    uv run zensical build

# Serve the documentation site locally
docs-serve: docs
    uv run zensical serve

# Deploy the documentation to GitHub Pages
docs-deploy: docs
    uv run zensical gh-deploy --force

# Set the version of all crates and python bindings
set-version VERSION:
    uv run python -c 'import re, sys; [open(f, "w").write(re.sub(r"(?m)^version\s*=\s*\".*\"", f"version = \"{sys.argv[1]}\"", content)) for f in ["Cargo.toml"] if (content := open(f).read())]' {{VERSION}}

# Run all Rust integration tests and Python pytest suite
test: sync
    cargo test --workspace
    uv run pytest tests/

# Run formatting and clippy linting checks
check:
    cargo fmt --all -- --check
    cargo clippy --workspace -- -D warnings

# Automatically format all Rust code
fmt:
    cargo fmt --all

# Run the Zani Rust engine benchmarks
bench:
    cargo bench --workspace

# Clean Cargo build artifacts and Python virtual environments
clean:
    cargo clean
    rm -rf .venv
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +
