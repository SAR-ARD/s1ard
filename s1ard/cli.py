import click


@click.command(name='s1ard',
               no_args_is_help=True,
               context_settings=dict(
                   ignore_unknown_options=True,
                   allow_extra_args=True, )
               )
@click.option('--config-file', '-c', required=False, type=click.Path(),
              help='Full path to an INI-style configuration text file.')
@click.option('--section', '-s', required=False, type=str, default='PROCESSING', show_default=True,
              help='Section of the configuration file to read processing related parameters from.')
@click.option('--debug', is_flag=True,
              help='Print debugging information for pyroSAR modules.')
@click.option('--version', is_flag=True,
              help='Print s1ard version information and exit. Overrides all other arguments.')
@click.pass_context
def cli(ctx, config_file, section, debug, version):
    """
    Central s1ard processing command.
    
    Additional options can be passed to override individual processing parameters
    in the configuration file. For example, to read all values from the configuration
    file except the acquisition mode and the annotation layers:
    
    s1ard -c config.ini --acq_mode IW --annotation dm,id
    
    The snap_gpt_args argument can be provided by using a dedicated argument separator and quotes:
    
    s1ard -c config.ini -- --snap_gpt_args "-J-Xmx100G -c 75G -q 30"
    
    \b
    The following defaults are set:
    (processing section)
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
    (metadata section)
    - access_url:        None
    - doi:               None
    - licence:           None
    - processing_center: None
    """
    import s1ard
    if version:
        print(s1ard.__version__)
    else:
        extra = {ctx.args[i][2:]: ctx.args[i + 1] for i in range(0, len(ctx.args), 2)}
        s1ard.process(config_file=config_file, section_name=section, debug=debug, **extra)
