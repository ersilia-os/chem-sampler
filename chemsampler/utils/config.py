import json


class Parameters(object):
    def __init__(self, params_json):
        with open(params_json, "r") as f:
            data = json.load(f)
        self.sampler_ids = data["samplers"]
        self.descriptor_ids = data["descriptors"]
        self.num_samples = int(data["num_samples"])
