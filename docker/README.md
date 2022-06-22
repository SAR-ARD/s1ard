# S1_NRB

# Docker Image

## Build

`docker build -t s1_nrb .`
## Run

`docker run -d -v [host_path_to_work_dir]:[container_path_to_work_dir] s1_nrb s1_nrb -c [path_to_config_file]`

- The [host_path_to_work_dir] path should contain all necessary inputs and output folders. 
  Paths written in the config file should resolve to the [container_path_to_work_dir].


