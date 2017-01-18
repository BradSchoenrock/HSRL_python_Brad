import re
import unicodedata
from enaml.validation import AbstractValidator
from enaml.validation import IntValidator, FloatValidator
from types import UnicodeType

class Validator(AbstractValidator):
    """Validator class inheriting AbstractValidator which verifies input values.
    
    The validator class is supplied directly to fields in the enaml file and they
    automatically execute it when a new value is enetered to verify whether or not
    it is valid input. It has been modified from the normal validator to change
    value types dynamically.
    """
    int_regex = re.compile('^[-]?\d+$')
    float_regex = re.compile('^[-]?\d+\.\d+$')
    float_okay_regex = re.compile('^[-]?\d*\.\d*$')
    str_regex = re.compile('^[a-zA-Z_]+$')
    okay_regex = re.compile('^$')
    input_type = "" 
        
    def set_type(self, input_type):
        """This method updates the value type being evaluated
        
        This method is called by the enaml file when a new attribute is selected
        and updates the value type accordingly which changes what type of regex
        is used in evaluation.
        """
        self.input_type = input_type
        
    def _validate(self, text):
        """This method validates the input and is called automatically by enaml."""
        
        #changes the text to unicode to ease evaluation.
        text = unicode(text)
        if self.input_type == "int":
            match = self.int_regex.match(text)
            if match:
                #print "good"
                return (self.ACCEPTABLE, text)
            
            match = self.okay_regex.match(text)
            if match:
                return (self.INTERMEDIATE, text)            
            else:
                #print "bad"
                return (self.INVALID, text)
            
        if self.input_type == "float":
                    match = self.float_regex.match(text)
                    if match:
                        #print "good"
                        return (self.ACCEPTABLE, text)
                    match = self.float_okay_regex.match(text)
                    if match:
                        return (self.INTERMEDIATE, text)                      
                    else:
                        #print "bad"
                        return (self.INVALID, text)
                    
        if self.input_type == "str" or self.input_type == "unicode":
                    match = self.str_regex.match(text)
                    if match:
                        #print "good"
                        return (self.ACCEPTABLE, text)
                    match = self.okay_regex.match(text)
                    if match:
                        return (self.INTERMEDIATE, text)                      
                    else:
                        #print "bad"
                        return (self.INVALID, text)
                    
        else:
            #print "bad"
            return (self.INVALID, text)            
    
    def validate(self,text):
        """This method is autocally called by the enaml file."""
        return self._validate(text)[0]
    
    def convert(self, text):
        """This method is autocally called by the enaml file."""
        return self._validate(text)[1]
    
    def format(self, text):
        """This method is autocally called by the enaml file."""
        return self._validate(text)[1]