FROM ubuntu:24.04 AS run-env

WORKDIR /root
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libgomp1 \
        libhts3t64 \
        libnuma1 \
        libboost-filesystem1.83.0 \
        libboost-graph1.83.0 \
        libboost-iostreams1.83.0 \
        libboost-program-options1.83.0 \
        libboost-regex1.83.0 \
        libboost-serialization1.83.0 \
        libboost-system1.83.0 \
        perl \
        python3 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

FROM run-env AS builder

ARG RABBITBIN_SOURCE_REVISION=unknown

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        autoconf \
        automake \
        git \
        libboost-all-dev \
        libdeflate-dev \
        libtool \
        cmake \
        libncurses-dev \
        pkg-config \
        zlib1g-dev \
        libhts-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY . /src/RabbitBin
RUN cd /src/RabbitBin && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DRABBITBIN_SOURCE_REVISION="${RABBITBIN_SOURCE_REVISION}" \
          -DRABBITBIN_NATIVE_ARCH=OFF .. && \
    make -j$(nproc) rabbitbin rabbit_depth rabbit_overlap && \
    make install && \
    rm -rf build

FROM run-env
ENV PATH=/usr/local/bin:$PATH
COPY --from=builder /usr/local /usr/local
CMD ["/usr/local/bin/run_rabbitbin.sh"]
