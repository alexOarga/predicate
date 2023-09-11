import sys
import cobra

def compute_name(name):
    return name + "[end]"

class NullDevice():
    def write(self, s):
        pass


def read_model(path):
    #original_stderr = sys.stderr  # keep a reference to STDERR
    #sys.stderr = NullDevice()  # redirect the real STDERR
    if path[-4:] == ".xml":
        model = cobra.io.read_sbml_model(path)
    elif path[-5:] == ".json":
        model = cobra.io.load_json_model(path)
    elif path[-4:] == ".yml":
        model = cobra.io.load_yaml_model(path)
    else:
        #sys.stderr = original_stderr  # turn STDERR back on
        raise RuntimeError("Model file must be either .xml .json .yml")
    #sys.stderr = original_stderr  # turn STDERR back on
    return model
