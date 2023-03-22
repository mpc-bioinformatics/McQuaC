# syntax=docker/dockerfile:1
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive
ENV PIPENV_YES=true
ENV PYTHONDONTWRITEBYTECODE 1

# Packages needed for MS-GF+ (java), TRFP (mono), installation stuff (curl, pipenv, pyenv) and libs needed for OpenMS 
RUN apt-get update && \
    apt-get install -y bash mono-complete openjdk-17-jre python3 python3-pip git curl make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
    libboost-date-time-dev libboost-iostreams-dev libboost-regex-dev libboost-math-dev \
    libboost-random-dev zlib1g libbz2-dev libsvm3 libxerces-c-dev libglpk-dev libqt5network5 libqt5opengl5 libqt5svg5 libqt5webkit5 libqt5core5a \
    libqt5sql5 && \
    apt-get clean
RUN curl https://pyenv.run | bash
RUN python3 -m pip install --user pipenv


# Export installed pipenv from python3
ENV PATH="${PATH}:/root/.local/bin" 

# Copy all the executables
COPY . /root/nextflowqc
WORKDIR /root/nextflowqc

# Install OpenMS and set up Pyenv
RUN chmod +x bin/* setup_pyenv_downlaod_extract_openms.sh
RUN bash setup_pyenv_downlaod_extract_openms.sh
# We can specify the needed nextflow version here!
RUN pipenv run pip install nextflow

# Export paths for executables from OpenMS
ENV PATH="${PATH}:/root/nextflowqc/bin/openms/usr/bin/" 

# Entrypoints and priotize pipenv executables
ENTRYPOINT PATH="$(pipenv --venv)/bin:${PATH}" bash
CMD PATH="$(pipenv --venv)/bin:${PATH}" bash
