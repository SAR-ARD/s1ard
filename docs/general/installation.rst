############
Installation
############

SNAP
----

S1_NRB requires ESA’s Sentinels Application Platform (SNAP) software to produce S1-NRB products. For download and installation
instructions please refer to the `official webpage <https://step.esa.int/main/download/snap-download/>`_.

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

    pip install git+https://github.com/SAR-ARD/S1_NRB.git@v0.4.1
