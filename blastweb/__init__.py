from flask import Flask
from .app import register_routes
import yaml
import os
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.middleware.proxy_fix import ProxyFix

def create_app(config_path=None, use_dispatcher=True):
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
    if use_dispatcher:
        prefix = config.get("url_prefix", "").rstrip("/")
        if prefix:
            app = DispatcherMiddleware(Flask("dummy"), {prefix: app})
    return app
