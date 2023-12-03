FROM --platform=amd64 mambaorg/micromamba:1.5.1-jammy
LABEL maintainer="dirk.winkelhardt@rub.de"

# Note: One could use a basic Python image to install the backend package
# but you have to install the mandatory dependencies via the OS package manager next to the conda environment
# making the package less maintainable.
# Another benefit is the use of a dedicated mamba user to adjust the UID/GID so the uploaded data is not owned by root.


ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000

USER root
RUN usermod "--login=mambauser" "--home=/home/mambauser" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=mambauser" \
        "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    # Update the expected value of MAMBA_USER for the
    # _entrypoint.sh consistency check.
    echo "mambauser" > "/etc/arg_mamba_user" && \
    :

# RUN apt-get update -y \
#     && apt-get install -y bash mono-complete openjdk-17-jre git curl make build-essential libssl-dev zlib1g-dev \
#     libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm unzip \
#     libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
#     libboost-date-time-dev libboost-iostreams-dev libboost-regex-dev libboost-math-dev \
#     libboost-random-dev zlib1g libbz2-dev libsvm3 libxerces-c-dev libglpk-dev libqt5network5 libqt5opengl5 libqt5svg5 libqt5webkit5 libqt5core5a \
#     libqt5sql5 \
#     && apt-get clean

RUN apt-get update -y \
    && apt-get install -y --no-upgrade git

WORKDIR /home/mambauser
# Copy backend and environment.yml
COPY --chown=mambauser:mambauser environment.yml .
COPY --chown=mambauser:mambauser requirements.txt .

USER mambauser
ENV HOME /home/mambauser
ENV PATH $PATH:$HOME/.local/bin
ENV ENV_NAME=nextqcflow


RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
