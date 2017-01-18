import json
import locate_file as lf
from collections import OrderedDict

def load_canvas_info(forModule=None):
	return json.load(file(lf.locate_file('artist_canvas.json',forModule=forModule)), object_pairs_hook=OrderedDict)
