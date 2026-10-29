from flask import Flask
from .app import register_routes, load_config
import os
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.middleware.proxy_fix import ProxyFix

def create_app(config_path=None, use_dispatcher=True):
    app = Flask(__name__)
    if config_path is None:
        config_path = os.environ.get("BLASTWEB_CONFIG", "blast.yaml")
    try:
        config = load_config(config_path)
    except (FileNotFoundError, ValueError) as e:
        raise RuntimeError(f"Invalid configuration: {e}")

    app.config.update(config)
    register_routes(app)
    if use_dispatcher:
        prefix = config.get("url_prefix", "")
        if prefix:
            app = DispatcherMiddleware(Flask("dummy"), {prefix.rstrip("/"): app})
    return app
