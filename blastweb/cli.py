import click
from blastweb.app import app

@click.command()
@click.option('--host', default='0.0.0.0', help='Host address')
@click.option('--port', default=5000, help='Port')
@click.option('--debug', is_flag=True, help='デバッグモードを有効にする')
def runserver(host, port, debug):
    """BLAST Web UI を起動します。"""
    app.run(host=host, port=port, debug=debug)

