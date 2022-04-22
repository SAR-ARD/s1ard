############
Installation
############

The S1_NRB package is not yet available via conda-forge or other common package distribution channels. In the meantime,
the following shall provide a convenient option provided that Anaconda or Conda has been installed:


1. Install and then activate the conda environment

```bash
conda env create --file https://raw.githubusercontent.com/SAR-ARD/S1_NRB/main/environment.yaml
conda activate nrb_env
```

2. Install the S1_NRB package into the environment

The package version can be changed as necessary (e.g., `@v0.2.0`).
See :ref:`changelog` for available versions.

```bash
pip install git+https://github.com/SAR-ARD/S1_NRB.git@v0.3.0
```
