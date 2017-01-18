import json
from collections import OrderedDict

missingnotset="missing isnt set blorp"

class json_config:
    """This class contains manages all attributes in the config file.
    
    This class is called when the configuration file is imported using the gui
    and creates a dictionary of all values in it. The dictionary contains two top
    level keys (more can be added) 'config' and 'display_defaults'. These two keys
    each contain subdictionaries which can be accessed by the methods of this class.
    All methods assume the display_defaults dictionary is being queried unless the key
    is specified as it is the most frequently accessed dictionary.
    """
    
    def __init__(self, cfg,default_key,allow_missing_values=False,allow_missing_enables=True):
        """init method which loads in the json file when class is called"""
        self.allow_missing_values=allow_missing_values
        self.allow_missing_enables=allow_missing_enables
        if default_key is None:
            raise RuntimeError("JSON Config default key cannot be None")
        try:
            if isinstance(cfg,basestring):
                self.confg_file = cfg
                f = open(cfg)
                self.obj = json.load(f, object_pairs_hook=OrderedDict)
                f.close()
            elif isinstance(cfg,dict):
                self.obj = cfg
                self.confg_file = None
            else:
                self.obj = cfg
                self.confg_file = cfg.config_file
        except IOError:
            raise #RuntimeError, "json_import : Could not read '%s' config file" %  cfg
        except ValueError as e:
            if e.message == "No JSON object could be decoded":#This is caught ONLY because the JSON parser doesn't give a useful error
                print 'Value error occurred reading JSON file',cfg,'. Possibly floating point value starts with decimal point, and needs a 0 prefix'
                raise ValueError('Value error occurred reading JSON file '+cfg+'. Possibly floating point value starts with decimal point, and needs a 0 prefix. '+repr(e))
            raise
        self.default_key=default_key
        if default_key not in self.obj:
            raise KeyError, "json_config can not find the required key '%s' in '%s'" % (default_key, cfg)
        self.usedAttrs={}

        # check if we're reading a old format configuration file
        #k = self.obj['display_defaults'].keys()[0]
        #if type(self.obj['display_defaults'][k]) == type([]):
        #    raise RuntimeError, 'Sorry - %s is an old-style configuration file' % cfg

    def copy(self):
        return json_config(self.obj,self.default_key,self.allow_missing_values,self.allow_missing_enables)

    def _setUsed(self,attr,key=None):
        if key is None:
           key=self.default_key
        if key not in self.usedAttrs:
            self.usedAttrs[key]=set()
        self.usedAttrs[key].add(attr)

    def json_representation(self):
        return self.obj
        
    def get_keys(self):
        """This method returns the top level keys when called.
        
        When called, this method will return a list of all the top level keys
        in the json file. This will normally be "config" and "display_defaults"
        but more can always be added in the future.
        This method is typically only called by other methods of this class.
        """
        
        return self.obj.keys()
    
    def get_labels(self,attr,key=None):
        """This method returns all labels for a given attribute.
        
        When called, this class will return all the labels associated with a
        given attribute, e.g., low, high, enable, min, max. If an incorrect 
        attribute is given, it will notify the user.
        """
        if key is None:
           key=self.default_key        
        try:
            ret=self.obj[key][attr].keys()
            self._setUsed(attr,key)
            return ret;
        except KeyError:
            if self.allow_missing_values:
                self._setUsed(attr,key)
                return tuple()
            raise KeyError("Could not find {0} in {1} ['{2}']".format(attr,key,self.confg_file))
        
    def get_attrs(self, key=None):
        """Returns all attributes for a given key when called.
        
        This method will return all attributes for a key, e.g. 'display_defaults'
        when called. This is primarily used for populating the list of attributes
        in the configuration tab.
        """
        if key is None:
           key=self.default_key        
        return self.obj[key].keys()
    
    def get_value(self,attr, label,key=None,require=False,missing=missingnotset):
        """Returns the value of a given label of an attribute.
        
        When called, this method returns the value of the label given
        for the given attribute. If the attribute or the label is invalid
        the user will be notified.
        """
        if key is None:
           key=self.default_key        
        try:
            ret= self.obj[key][attr][label]
            self._setUsed(attr,key)
            return ret;
        except KeyError:
            if (missing is not missingnotset or self.allow_missing_values) and (not require or missing is not missingnotset):
                self._setUsed(attr,key)
                if missing is missingnotset:
                    return None
                return missing
            raise KeyError("Could not find {0} or {1} in {2} ['{3}']".format(attr,label,key,self.confg_file))

            
    def set_value(self,attr, label, value,key=None):
        """This method changes the value for a given label in an attribute.
        
        When called, this method takes the value that is passed to it for a 
        given label of an attribute and sets it as the new value. This is 
        primarily called when saving the configuration file
        """
        if key is None:
           key=self.default_key        
        self._setUsed(attr,key)
        if attr not in self.obj[key]:
            self.obj[key][attr]=OrderedDict()
        self.obj[key][attr][label] = value
            
    def get_all_data(self):
        """Returns all data for a given key in the form of nested lists.
        
        When called, this method will get all data available for every attribute
        for the given key and return it as a nested list. This is primarily used
        by configModel for managing the data for the gui.
        """
        
        data = []
        for attr in self.get_attrs():
            temp = []
            for label in self.get_labels(attr):
                temp.append(self.get_value(attr,label))
            data.append(temp)
        return data
    
    def get_size(self, attr):
        """Returns a list of the x and y values of an attribute
        which contains those labels
        
        This method returns the x and y values in a list of an attribute that
        has those labels. This is considered a special method of this class
        as it is primarily used by image_size as to not complicate the code
        in display_utilities.py
        """
        x = self.get_value(attr, 'X',require=True)
        y = self.get_value(attr, 'Y',require=True)
        return (x,y)
        
    def enabled(self, attr, key=None, return_if_missing=missingnotset):
        """Returns 1 if the attribute exists and enable is set to 1, else 0
        
        When called, this method checks the following:
          has a given key (optional) and attr, or will return False
          attr has 'enable', or will return True
          returns truth value of value at 'enable'
        """
        if key is None:
            key=self.default_key
        if return_if_missing is not missingnotset:
            assert(isinstance(return_if_missing,bool))
        if attr not in self.get_attrs(key) and self.allow_missing_enables:#entire attr is missing
            self._setUsed(attr,key)
            if return_if_missing is missingnotset:
                return False
            return return_if_missing
        try:
            if self.get_value(attr, 'enable', key):#has truth value of true. no numbers or objects will be returned by this function
                return True
            return False
        except KeyError: #just enable is missing, attr exists
            if self.allow_missing_enables:
                self._setUsed(attr,key)
                if return_if_missing is missingnotset:
                    return True
                return return_if_missing
            raise
   

    def items(self):
        """This method emulates what .items() returned with the old configuration.
        
        This method will format the data for a given key into a nested lists where each
        entry is formated as such: [attr [list of values]]. This exists to emulate how
        calling the items() function worked on the old dictionaries as to not complicate
        code in rti.py. See rti.py ~line 800.
        """
        
        attr = self.get_attrs()
        data = self.get_all_data()
        items = []
        for i in range(len(attr)):
            items.append([attr[i], data[i]])
        return items
    
    def get_all_values(self,attr,key=None):
        """Returns all given values for the supplied attribute.
        
        When called, this method returns all the values of a given attribute
        as a list. This is a special case method used by display_utilities.py for
        getting the values of the attribute counts vs time as it can be expanded
        indefinitely with the number of labels associated with it and is usually
        checks the 1st indice to the end indice [1:].
        """
        if key is None:
            key=self.default_key
        ret=self.obj[key][attr].values()
        self._setUsed(attr,key)
        return ret
        
    def get_dict(self, key=None):
        """ returns a dictionary included in attr"""
        return self.obj[key or self.default_key]   

    def get_filtereddict(self, key=None):
        """ returns a dictionary included in attr"""
        if key is None:
            key=self.default_key
        used=self.obj[key].copy()
        unused=OrderedDict()
        if key in self.usedAttrs:
            tmp=used
            used=OrderedDict()
            for k in tmp.keys():
                if k in self.usedAttrs[key]:
                    used[k]=tmp[k]
                else:
                    unused[k]=tmp[k]
        return used,unused

import lg_base.core.locate_file as lf

def get_display_defaults(display_defaults_file):
    """get display defaults from display_default.json.
       image selections"""
    display_defaults = json_config(lf.locate_file(display_defaults_file),'display_defaults',allow_missing_values=True)
    config = display_defaults
    return (display_defaults, config)



def main():
    import sys
    v=json_config(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'unknown')

if __name__ == '__main__':
    main()
