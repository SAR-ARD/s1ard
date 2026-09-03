Installation
============

SNAP
----

s1ard offers an interface to plug in various SAR processing software solutions.
The plugin for ESA's Sentinels Application Platform (SNAP) is directly included in s1ard.
Downloaders for different operating systems can be obtained from the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

The following code can be used to replicate the software installation on a Linux OS:

::

    VERSION=13
    TARGET=~/SNAP"$VERSION"

    INSTALLER=esa-snap_sentinel_linux-"$VERSION".0.0.sh
    wget https://download.esa.int/step/snap/"$VERSION".0/installers/"$INSTALLER"
    bash $INSTALLER -q -dir $TARGET
    $TARGET/bin/snap --nosplash --nogui --modules --update-all

    # add SNAP location to the PATH environment variable in the .bashrc file
    echo PATH=$PATH:$TARGET/bin >> ~/.bashrc

See also the web page on how to `update SNAP from the command line <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line>`_.

Alternatively, updates for individual modules and versions can be downloaded in the `SNAP Update Center <https://step.esa.int/updatecenter/>`_.

s1ard
------

The s1ard package is already available via PyPI.
However, the installation via mamba is currently preferred, in particular due to the GDAL dependency.
The following shall provide a convenient installation option provided that Anaconda or Miniconda has been installed.

Latest State on Github
++++++++++++++++++++++

1. Create and then activate the conda environment

::

    mamba env create --file https://raw.githubusercontent.com/SAR-ARD/s1ard/main/environment.yml
    mamba activate s1ard

2. Install the s1ard package into the environment

::

    pip install git+https://github.com/SAR-ARD/s1ard.git

Specific Version
++++++++++++++++

The package version can be changed as necessary. See the `Tags <https://github.com/SAR-ARD/s1ard/tags>`_ section of the
repository for available versions.

::

    # in older versions the file was called environment.yaml
    mamba env create --file https://raw.githubusercontent.com/SAR-ARD/s1ard/v1.0.0/environment.yaml
    mamba activate s1ard
    pip install git+https://github.com/SAR-ARD/s1ard.git@v1.0.0

Docker
------

Both SNAP and s1ard can also be installed into a docker container using the Dockerfile that is provided with the package.
