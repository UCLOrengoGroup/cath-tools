# This Dockerfile provides an illustrative example of how to build cath-tools

FROM ubuntu:20.04

# DEBIAN_FRONTEND=noninteractive prevents the tzdata package hanging on an interactive user prompt
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    cmake \
    g++-10 \
    git \
    libgsl-dev \
    ninja-build \
    python3-pip
RUN /usr/bin/pip3 install --upgrade conan

RUN conan profile new default --detect && conan profile update settings.compiler.libcxx=libstdc++11 default

WORKDIR /cath-data
ARG BUILD_EXTRA_CATH_TESTS=ON
ARG BUILD_EXTRA_CATH_TOOLS=ON

ADD . /cath-data
RUN conan install --build missing --install-folder docker-build-dir . --profile .github/conan-profiles/ubuntu-20.04.gcc.release.conan-profile
RUN cmake -GNinja -B docker-build-dir -S . -D CMAKE_TOOLCHAIN_FILE:FILEPATH=.github/cmake-toolchain-files/ubuntu-20.04.gcc.release.cmake -D CMAKE_MODULE_PATH:PATH=/cath-data/docker-build-dir
RUN ninja -C docker-build-dir -k 0
RUN ( cd docker-build-dir && ctest -j 6 --output-on-failure )
