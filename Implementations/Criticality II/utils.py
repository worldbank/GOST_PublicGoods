import os,json

def load_config():

    """Read config.json

    """

    config_path = os.path.join(os.path.dirname(__file__), '..', 'config.json')

    with open(config_path, 'r') as config_fh:

        config = json.load(config_fh)

    return config
