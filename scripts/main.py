import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, ".."))
from chemsampler.runner import Runner

config_file = os.path.abspath(os.path.join(root, "..", "default", "params.json"))

#Runner(config_file=config_file).run()
Runner(config_file=config_file).add_properties()
