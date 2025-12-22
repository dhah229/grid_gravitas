# -----------------
# Stage 1: Builder
# -----------------
FROM rust:slim-bookworm as builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    pkg-config \
    libclang-dev \
    clang \
    cmake \
    sqlite3 \
    libsqlite3-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/grid_gravitas
COPY . .

RUN cargo build --release

# -----------------
# Stage 2: Runtime
# -----------------
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y \
    libgdal32 \
    libgeos-c1v5 \
    libproj25 \
    libnetcdf19 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/src/grid_gravitas/target/release/grid_gravitas /usr/local/bin/grid_gravitas

WORKDIR /data

ENTRYPOINT ["grid_gravitas"]

CMD ["--help"]
