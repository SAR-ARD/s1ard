# S1-NRB prototype processor

## Installation

- Install and then activate conda environment: 
```bash
conda env create --file environment.yaml
conda activate nrb_env
```

- Install this package into the environment:
```bash
pip install .
```

## Usage

- Create a config file (see `config.ini` for an example) and store it in your project directory
- Print a help message using:
```bash
s1_nrb --help
```

- Start the processor using parameters defined in a config file:
```bash
s1_nrb -c /path/to/your/config.ini
```

A description of individual steps can be found [here](https://github.com/SAR-ARD/S1_NRB/blob/main/docs/S1_NRB.rst).
