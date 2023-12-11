from chemsampler.utils.config import ConfigRun

config_file = "default/params.json"
rc = ConfigRun(config_file)
rc.create_output_files()
params = rc.read_config_file
print(params)


