Installation
============

SNAP
----

S1_NRB requires ESA’s Sentinels Application Platform (SNAP) software to produce S1-NRB products. It has been developed based on SNAP version 8.
Downloaders for different operating systems can be obtained from the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

The following code can be used to replicate the software installation on a Linux OS:

::

    VERSION=8
    TARGET=~/SNAP"$VERSION"

    INSTALLER=esa-snap_sentinel_unix_"$VERSION"_0.sh
    wget https://download.esa.int/step/snap/"$VERSION".0/installers/"$INSTALLER"
    bash $INSTALLER -q -dir $TARGET
    $TARGET/bin/snap --nosplash --nogui --modules --update-all

See also the web page on how to `update SNAP from the command line <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line>`_.

Alternatively, updates for individual modules and versions can be downloaded in the `SNAP Update Center <https://step.esa.int/updatecenter/>`_.

S1_NRB
------

The S1_NRB package is not yet available via conda-forge or other common package distribution channels. In the meantime,
the following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed:

1. Create and then activate the conda environment

::

    conda env create --file https://raw.githubusercontent.com/SAR-ARD/S1_NRB/main/environment.yaml
    conda activate nrb_env

2. Install the S1_NRB package into the environment

The package version can be changed as necessary. See the `Tags <https://github.com/SAR-ARD/S1_NRB/tags>`_ section of the
repository for available versions.

::

    pip install git+https://github.com/SAR-ARD/S1_NRB.git@v1.0.0
