## These commands assume that you have the necessary dependencies installed and configured for each language.

# Run the astroz (Zig) SGP4 benchmark
bench-astroz-zig:
    zig build bench -Doptimize=ReleaseFast

# Run the Rust sgp4 benchmark
bench-rust:
    rust-script benchmarks/rust_bench.rs

# Run the Python sgp4 benchmark
bench-python-sgp4:
    uv run benchmarks/python_sgp4_bench.py

# Run the Python astroz benchmark
bench-astroz-python:
    uv run benchmarks/python_astroz_bench.py

# Run the JAX CPU SGP4 benchmark
bench-jax-cpu:
    uv run benchmarks/jax_cpu_bench.py

# Run the JAX GPU SGP4 benchmark
bench-jax-gpu:
    uv run benchmarks/jax_gpu_bench.py
