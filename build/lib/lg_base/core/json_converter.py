#!/usr/bin/env python
from json_config import json_config
import json
import os

def convert_to_new(master, original ):    
    """
    convert old style plot selection files to new format    
    """
    json_file = json_config(master)
    f = open(original)
    old_defaults = json.load(f)
    f.close()
    if not old_defaults.has_key('display_defaults'):
        print 'skipping  %s - not a display configuration' % original
        return
    # iterate through all plot_selections, like 'backscat_image'
    for attr in old_defaults['display_defaults'].keys():
        # skip old comments
        if not attr.startswith('#'):
            values = old_defaults['display_defaults'][attr]
            i = 0
            try:
                for label in json_file.get_labels(attr):
                    # if the original file didn't contain all the values,
                    # use the values from the master, instead.
                    if i == len(values):
                        print 'Note: %s had missing values for %s' % (original, attr)
                        break
                    val = values[i]
                    json_file.set_value(attr,label,val)
                    i = i+ 1
            except KeyError:
                print 'skipping attribute %s' % attr

    backup = original + '.old'
    os.rename(original, backup)
    data = json.dumps(json_file.obj, indent=4)
    f = open(original,'w')
    f.write(data)   
    f.close()

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print 'usage: %s master_configuration old_style_json' % sys.argv[0]
        sys.exit(1)

    for i in range(2,len(sys.argv)):
        filename = sys.argv[i]
        print 'converting %s ' % filename
        convert_to_new(sys.argv[1], filename)
