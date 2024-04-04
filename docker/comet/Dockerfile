FROM --platform=amd64 ubuntu:22.04

ARG COMET_VERSION=v2022.01.2

ENV DEBIAN_FRONTEND noninteractive
ENV COMET_VERSION=${COMET_VERSION}

COPY docker/comet/entrypoint.sh /usr/local/bin

RUN apt-get update -y \
    && apt-get install -y --no-upgrade curl

RUN curl -L -o /usr/local/bin/comet https://github.com/UWPR/Comet/releases/download/${COMET_VERSION}/comet.linux.exe \
    && chmod +x /usr/local/bin/comet

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]