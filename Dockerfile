FROM mambaorg/micromamba:1.4.3-focal
# https://hub.docker.com/r/mambaorg/micromamba

USER root
COPY . /src/

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends --no-install-suggests \
    build-essential \
    fontconfig \
    fonts-dejavu \
    git \
    libgfortran5 \
    wget \
    zip \
    && apt-get autoremove -y \
    && apt-get purge -y --auto-remove \
    && rm -rf /var/lib/apt/lists/*

RUN bash /src/docker/esa-snap.sh
ENV PATH="/usr/local/snap/bin:${PATH}"

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -f /src/environment.yaml \
    && pip install /src/ \
    && micromamba clean -ay \
    && rm -rf /src/