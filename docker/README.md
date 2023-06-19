# Docker

## Build

After cloning the repository, you can build the Dockerfile with the following command. You can replace the `s1_nrb` tag with any other tag you want.  

    docker build -t s1_nrb .

Helpful `docker build` [options](https://docs.docker.com/engine/reference/commandline/build/#options):
- `--no-cache` builds the image without using cached layers. Might be useful if you encounter problems during the build process.
- `--platform` builds the image for a specific platform. E.g. `--platform linux/amd64` builds the image for the amd64 architecture.

## Run

    docker run -v [host_path_to_work_dir]:[container_path_to_work_dir] s1_nrb s1_nrb -c [path_to_config_file]

- `[host_path_to_work_dir]` should contain all necessary input data on the host system.
- `[container_path_to_work_dir]` should be the path to the work directory in the container. E.g. `/data`
- First `s1_nrb` is the given tag of the image (see comment in the build section), while the second `s1_nrb` is the name of the executable.
- `[path_to_config_file]` should be relative to the `[container_path_to_work_dir]` and paths written in the config file should resolve to the `[container_path_to_work_dir]`.

Example:

    docker run -v /home/user/work_dir:/data s1_nrb s1_nrb -c /data/config.ini
