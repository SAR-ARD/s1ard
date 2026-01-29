# s1ard

# Docker Image

## Build

From the s1ard directory run:
`bash docker/build_image.sh`

Possible options are:
- `-s` which installs SNAP inside the image
- `-n` which disables the Docker build cache

Using both the command looks like:
`bash docker/build_image.sh -s -n`

## Run
`docker run -d -v [host_path_to_work_dir]:[container_path_to_work_dir] s1ard-[base/snap]:[version] s1ard -c [path_to_config_file]`

- The [host_path_to_work_dir] path should contain all necessary inputs and output folders.
  Paths written in the config file should resolve to the [container_path_to_work_dir].
