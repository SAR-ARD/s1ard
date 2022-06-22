FROM conda/miniconda3

USER root

RUN apt-get update && apt-get install -y \
    software-properties-common
RUN apt-get install -y git python3-pip wget libpq-dev


WORKDIR /tmp/
RUN wget https://download.esa.int/step/snap/8.0/installers/esa-snap_sentinel_unix_8_0.sh
COPY docker/esa-snap.varfile /tmp/esa-snap.varfile
RUN chmod +x esa-snap_sentinel_unix_8_0.sh

## install and update snap
RUN /tmp/esa-snap_sentinel_unix_8_0.sh -q /tmp/varfile esa-snap.varfile
RUN apt install -y fonts-dejavu fontconfig
COPY docker/update_snap.sh /tmp/update_snap.sh
RUN chmod +x update_snap.sh
RUN /tmp/update_snap.sh

# install NRB
SHELL [ "/bin/bash", "--login", "-c" ]

COPY environment.yaml environment.yaml
RUN conda env create --force  --file environment.yaml

RUN echo "conda init bash" >> ~/.bashrc
RUN source ~/.bashrc
RUN  echo "conda activate nrb_env" >> ~/.bashrc

WORKDIR /app/
COPY . /app/
RUN source ~/.bashrc \
 && python setup.py install


COPY docker/entrypoint.sh entrypoint.sh
RUN chmod +x entrypoint.sh
ENTRYPOINT ["/app/entrypoint.sh"]





