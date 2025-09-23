import click


@click.group(name='s1rb',
             no_args_is_help=True,
             invoke_without_command=True)
@click.option('--version', is_flag=True,
              help='Print s1ard version information and exit. Overrides all other arguments.')
def cli(version):
    if version:
        import s1ard
        print(s1ard.__version__)


@cli.command(name='init',
             no_args_is_help=True,
             context_settings=dict(
                 ignore_unknown_options=True,
                 allow_extra_args=True, )
             )
@click.option('--config-file', '-c', required=True, type=click.Path(),
              help="Full path to an INI-style target configuration text file.")
@click.option('--overwrite', '-o', is_flag=True, default=False,
              help='Overwrite an existing file?')
@click.option('--config-source', '-s', required=False, type=click.Path(),
              help="Full path to an INI-style source configuration text file. "
                   "If not defined, configuration will be read from the package's default file.")
@click.pass_context
def init(ctx, config_file, overwrite=False, config_source=None):
    """
    Initialize configuration for s1ard radar backscatter processing.

    This creates a processor configuration file in a user-defined location.
    The package's default file or a user-defined source file may serve as base.
    Additional options can be passed to override individual processing parameters
    in the config file. For example, to read all values from the default
    file except the acquisition mode and the annotation layers:

    s1rb init -c config.ini --acq_mode IW --annotation dm,id
    
    A previously created config file can be used instead of the default file:
    
    s1rb init -c config.ini -s config_source.ini --acq_mode IW --annotation dm,id

    SAR processor command line arguments (like SNAP gpt_args) can be provided
    by using a dedicated argument separator and quotes:

    s1rb init -c config.ini -- --gpt_args "-J-Xmx100G -c 75G -q 30"

    \b
    The following defaults are set:
    PROCESSING
    - annotation:        dm,ei,id,lc,li,np,ratio
    - dem_type:          Copernicus 30m Global DEM
    - date_strict:       True
    - etad:              False
    - etad_dir:          None
    - gdal_threads:      4
    - measurement:       gamma
    - ard_dir:           ARD
    - sar_dir:           SAR
    - tmp_dir:           TMP
    - wbm_dir:           WBM
    - logfile:           None
    METADATA
    - access_url:        None
    - doi:               None
    - licence:           None
    - processing_center: None
    """
    from s1ard.config import init as init_config
    extra = {ctx.args[i][2:]: ctx.args[i + 1]
             for i in range(0, len(ctx.args), 2)}
    init_config(target=config_file, source=config_source,
                overwrite=overwrite, **extra)


@cli.command(name='process',
             no_args_is_help=True,
             context_settings=dict(
                 ignore_unknown_options=True,
                 allow_extra_args=True, )
             )
@click.option('--config-file', '-c', required=False, type=click.Path(),
              help="Full path to an INI-style configuration text file. "
                   "If not defined, the package's default file will be used.")
@click.option('--debug', is_flag=True,
              help='Print debugging information for pyroSAR modules.')
@click.pass_context
def process(ctx, config_file, debug):
    """
    Central s1ard radar backscatter (rb) processing command.
    
    A config file can be defined from which all configuration is read.
    If not defined, configuration will be read from the package's default file.
    Additional options can be passed to override individual processing parameters
    in the configuration file. For example, to read all values from the configuration
    file except the acquisition mode and the annotation layers:
    
    s1rb process -c config.ini --acq_mode IW --annotation dm,id
    
    The snap_gpt_args argument can be provided by using a dedicated argument separator and quotes:
    
    s1rb process -c config.ini -- --snap_gpt_args "-J-Xmx100G -c 75G -q 30"
    
    \b
    The following defaults are set:
    PROCESSING
    - annotation:        dm,ei,id,lc,li,np,ratio
    - dem_type:          Copernicus 30m Global DEM
    - date_strict:       True
    - etad:              False
    - etad_dir:          None
    - gdal_threads:      4
    - snap_gpt_args:     None
    - measurement:       gamma
    - ard_dir:           ARD
    - sar_dir:           SAR
    - tmp_dir:           TMP
    - wbm_dir:           WBM
    - logfile:           None
    METADATA
    - access_url:        None
    - doi:               None
    - licence:           None
    - processing_center: None
    """
    from s1ard import process
    extra = {ctx.args[i][2:]: ctx.args[i + 1] for i in range(0, len(ctx.args), 2)}
    process(config_file=config_file, debug=debug, **extra)
