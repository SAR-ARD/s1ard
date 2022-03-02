import click
import S1_NRB.processor as process


@click.command()
@click.option('--config-file', '-c', required=True, type=click.Path(),
              help='Full path to an INI-style configuration text file.')
@click.option('--section', '-s', required=False, type=str, default='GENERAL', show_default=True,
              help='Section of the configuration file to read parameters from.')
@click.option('--debug', is_flag=True,
              help='Print debugging information for pyroSAR modules.')
def cli(config_file, section, debug):
    process.main(config_file=config_file, section_name=section, debug=debug)
