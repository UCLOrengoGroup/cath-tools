# This Dockerfile provides an illustrative example of how to build cath-tools

FROM ubuntu:20.04

ENV TZ="UTC"

RUN apt-get update && apt-get install -y apt-transport-https ca-certificates gnupg software-properties-common wget
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ focal main'

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