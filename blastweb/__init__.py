from flask import Flask
from .app import register_routes
import yaml

def create_app(config_path=None):
    app = Flask(__name__)
    if config_path is None:
        config_path = os.environ.get("BLASTWEB_CONFIG", "blast.yaml")
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        raise RuntimeError(f"Configuration file not found: {config_path}")

    app.config.update(config)
    register_routes(app)

    return app
