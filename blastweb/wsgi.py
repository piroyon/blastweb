from blastweb import create_app
import os

config_path = os.environ.get("BLASTWEB_CONFIG", "blast.yaml")
app = create_app(config_path)

