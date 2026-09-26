# Benchmark recipes. These assume each language's toolchain is installed
# (zig, rust-script, uv).

# Run the astroz Zig SGP4 benchmark
bench-astroz-zig:
    zig build bench -Doptimize=ReleaseFast

# Run the Rust sgp4 crate benchmark
bench-rust:
    rust-script benchmarks/rust_bench.rs

# Run the python-sgp4 benchmark
bench-python-sgp4:
    uv run benchmarks/python_sgp4_bench.py

# Run the astroz Python bindings benchmark
bench-astroz-python:
    uv run benchmarks/python_astroz_bench.py

# Run the astrojax benchmark on CPU
bench-jax-cpu:
    uv run benchmarks/jax_cpu_bench.py

# Run the astrojax benchmark on GPU
bench-jax-gpu:
    uv run benchmarks/jax_gpu_bench.py
