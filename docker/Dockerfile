# syntax=docker/dockerfile:1
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive
ENV PIPENV_YES=true
ENV PYTHONDONTWRITEBYTECODE 1

# Packages needed for MS-GF+ (java), TRFP (mono), installation stuff (curl, pipenv, pyenv) and libs needed for OpenMS 
RUN apt-get update && \
    apt-get install -y bash mono-complete openjdk-17-jre python3 python3-pip git curl make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm python-is-python3 unzip \
    libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
    libboost-date-time-dev libboost-iostreams-dev libboost-regex-dev libboost-math-dev \
    libboost-random-dev zlib1g libbz2-dev libsvm3 libxerces-c-dev libglpk-dev libqt5network5 libqt5opengl5 libqt5svg5 libqt5webkit5 libqt5core5a \
    libqt5sql5 && \
    apt-get clean

# Copy all the executables
COPY . /root/nextflowqc
WORKDIR /root/nextflowqc

# Install dependencies and set up Pyenv
RUN chmod +x setup_pyenv_downlaod_extract_openms.sh
RUN bash setup_pyenv_downlaod_extract_openms.sh
RUN chmod +x bin/*

# We can specify the needed nextflow version here!
RUN python -m pip install nextflow

# Export paths for executables from OpenMS and other bin executables
ENV PATH="${PATH}:/root/nextflowqc/bin:/root/nextflowqc/bin/openms/usr/bin/" 

# Entrypoints and priotize pipenv executables
ENTRYPOINT bash
CMD bash
