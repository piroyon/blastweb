import os
import shutil
import click
from flask import Flask
from blastweb import create_app
import importlib.resources as pkg_resources

@click.group()
def cli():
    """BLAST Web UI Command Line Interface."""
    pass

@cli.command()
@click.option('--host', default='0.0.0.0', help='Host to bind')
@click.option('--port', default=5001, help='Port to bind')
@click.option('--debug', is_flag=True, help='Enable debug mode')
@click.option('--config', default='blast.yaml', help='Path to blast.yaml')
def runserver(host, port, debug, config):
    """Run the development server."""
    app = create_app(config)
    app.run(host=host, port=port, debug=debug)

@cli.command()
def init():
    """Create a blast.yaml from example in current directory."""
    try:
        example_path = pkg_resources.files("blastweb").joinpath("blast.yaml.example")
    except AttributeError:
        # Python <3.9 fallback
        import pkg_resources as legacy_pkg_resources
        example_path = legacy_pkg_resources.resource_filename("blastweb", "blast.yaml.example")

    dest_path = os.path.join(os.getcwd(), "blast.yaml")
    if os.path.exists(dest_path):
        click.echo("blast.yaml already exists in this directory.")
    else:
        shutil.copy(example_path, dest_path)
        click.echo("blast.yaml created.")


if __name__ == "__main__":
    cli()

