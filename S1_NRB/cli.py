import click


@click.command(name='s1_nrb',
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
              help='Print S1_NRB version information and exit. Overrides all other arguments.')
@click.pass_context
def cli(ctx, config_file, section, debug, version):
    """
    Central S1_NRB processing command.
    
    Additional options can be passed to override individual processing parameters
    in the configuration file. For example, to read all values from the configuration
    file except the acquisition mode and the annotation layers:
    
    s1_nrb -c config.ini --acq_mode IW --annotation dm,id
    """
    import S1_NRB
    if version:
        print(S1_NRB.__version__)
    else:
        extra = {ctx.args[i][2:]: ctx.args[i + 1] for i in range(0, len(ctx.args), 2)}
        S1_NRB.process(config_file=config_file, section_name=section, debug=debug, **extra)
