Installation
============

SNAP
----

s1ard requires ESA's Sentinels Application Platform (SNAP) software for SAR data processing.
Version 1.0.0 has been developed based on SNAP 8.
SNAP 9 is supported since version 1.0.2, SNAP 10 since version 1.7.0.
Downloaders for different operating systems can be obtained from the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

The following code can be used to replicate the software installation on a Linux OS:

::

    VERSION=10
    TARGET=~/SNAP"$VERSION"

    INSTALLER=esa-snap_sentinel_linux-"$VERSION".0.0.sh
    wget https://download.esa.int/step/snap/"$VERSION"_0/installers/"$INSTALLER"
    bash $INSTALLER -q -dir $TARGET
    $TARGET/bin/snap --nosplash --nogui --modules --update-all

    # add SNAP location to the PATH environment variable in the .bashrc file
    echo PATH=$PATH:$TARGET/bin >> ~/.bashrc

See also the web page on how to `update SNAP from the command line <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line>`_.

Alternatively, updates for individual modules and versions can be downloaded in the `SNAP Update Center <https://step.esa.int/updatecenter/>`_.
The latest bundle that was used during release of version 1.0.0 is available here: https://step.esa.int/updatecenter/8.0_20220323-143356/.

s1ard
------

The s1ard package is not yet available via conda-forge or other common package distribution channels. For now,
the following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed.

Latest State on Github
++++++++++++++++++++++

1. Create and then activate the conda environment

::

    conda env create --file https://raw.githubusercontent.com/SAR-ARD/s1ard/main/environment.yaml
    conda activate s1ard

2. Install the s1ard package into the environment

::

    pip install git+https://github.com/SAR-ARD/s1ard.git

Specific Version
++++++++++++++++

The package version can be changed as necessary. See the `Tags <https://github.com/SAR-ARD/s1ard/tags>`_ section of the
repository for available versions.

::

    conda env create --file https://raw.githubusercontent.com/SAR-ARD/s1ard/v1.0.0/environment.yaml
    conda activate s1ard
    pip install git+https://github.com/SAR-ARD/s1ard.git@v1.0.0

Docker
------

Both SNAP and s1ard can also be installed into a docker container using the Dockerfile that is provided with the package.
