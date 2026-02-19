### ----------------------------------------------- ###
### --             Base conda image              -- ###
### ----------------------------------------------- ###
FROM condaforge/miniforge3:25.11.0-0 AS conda-base

# Installing prerequisites
RUN apt update && \
    apt upgrade -y && \
    apt install -y \
        git \
        wget \
        libxext-dev \
        libxrender1 \
        libxtst6 \
        libxi6 \
    && apt clean

### ----------------------------------------------- ###
### --             Base s1ard image              -- ###
### ----------------------------------------------- ###
FROM conda-base AS s1ard-base

RUN mkdir /app
COPY . /app
WORKDIR /app

# Create conda environment and install s1ard
RUN conda env create --yes --file /app/environment.yaml
RUN conda run -n s1ard python -m pip install .

# Set up general environment
ENV CONDA_DEFAULT_ENV=s1ard
ENV PATH=/opt/conda/envs/s1ard/bin:$PATH

### ---------------------------------------------------- ###
### --             s1ard with SNAP image              -- ###
### ---------------------------------------------------- ###
FROM s1ard-base AS s1ard-snap

# Set desired SNAP version
ENV SNAP_VERSION=13.0
ENV TARGET_DIR=/usr/local/bin/SNAP${SNAP_VERSION}

# Download the corresponding SNAP installer
WORKDIR /tmp
RUN wget https://download.esa.int/step/snap/${SNAP_VERSION}/installers/esa-snap_sentinel_linux-${SNAP_VERSION}.0.sh
RUN cp /app/docker/esa-snap.varfile /tmp/esa-snap.varfile
RUN chmod +x esa-snap_sentinel_linux-${SNAP_VERSION}.0.sh

# Install SNAP
RUN /tmp/esa-snap_sentinel_linux-${SNAP_VERSION}.0.sh -q /tmp/varfile esa-snap.varfile

# Update SNAP (see https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line)
SHELL ["/bin/bash", "-c"]
RUN ${TARGET_DIR}/bin/snap --nosplash --nogui --modules --update-all 2>&1 | \
    while read -r line; do \
        echo "$line"; \
        if [[ "$line" == "updates=0" ]]; then \
            sleep 2; \
            pkill -TERM -f "snap/jre/bin/java"; \
        fi; \
    done

WORKDIR /app
