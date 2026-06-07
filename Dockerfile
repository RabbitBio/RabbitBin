FROM ubuntu:24.04 AS run-env

WORKDIR /root
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends libgomp1 python3 libhts3 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

FROM run-env AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libboost-all-dev \
        cmake \
        libncurses5-dev \
        zlib1g-dev \
        libhts-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY . /src/RabbitBin
RUN cd /src/RabbitBin && \
    mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) rabbitbin && \
    make install && \
    rm -rf build

FROM run-env
ENV PATH=/usr/local/bin:$PATH
COPY --from=builder /usr/local /usr/local
CMD ["/usr/local/bin/run_rabbitbin.sh"]
