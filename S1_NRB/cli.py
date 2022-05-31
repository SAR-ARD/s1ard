import click
import S1_NRB


@click.command()
@click.option('--config-file', '-c', required=False, type=click.Path(),
              help='Full path to an INI-style configuration text file.')
@click.option('--section', '-s', required=False, type=str, default='GENERAL', show_default=True,
              help='Section of the configuration file to read parameters from.')
@click.option('--debug', is_flag=True,
              help='Print debugging information for pyroSAR modules.')
@click.option('--version', is_flag=True,
              help='Print S1_NRB version information. Overrides all other arguments.')
def cli(config_file, section, debug, version):
    if version:
        print(S1_NRB.__version__)
    else:
        S1_NRB.process(config_file=config_file, section_name=section, debug=debug)
