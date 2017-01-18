import re
import os
import json
from enaml.item_models.abstract_item_model import AbstractItemModel
from parsers import parse_json
from enaml.item_models.standard_models import ListModel
from traits.api import HasTraits, Str, Any, Bool, Int, Unicode, List, on_trait_change
from json_config import json_config
        
class FieldModel(HasTraits):
    """This class is called by each configForm present in the enaml file.
    
    This class contains the methods needed by the configForm class in enaml
    to update accordingly to the data type which they represent at the time.
    """
    #These assignments needed by enaml.
    label = Unicode
    value = Any
    visible = Bool
    value_type = Str
    checked = Bool
    
    def __init__(self):
        self.visible = False
        self.checked = False
        
    def set_label(self,label):
        """Allows the label of the form to be set"""
        
        self.label = label
        
    def set_value(self,value):
        """Allows the value of the form to be set"""
        
        self.value = value
        
    def set_visible(self,visible):
        """Allows the form to be set to visibile or invisible"""
        self.visible = visible
        
    def set_value_type(self,value_type):
        """Allows the value_type of the form to be set"""
        self.value_type = value_type
        
    def set_checked(self,checked):
        """Allows whether or not the form is checked to be set"""
        self.checked = checked
        
    def check_bool(self):
        """Checks to see if the label is "enable" and if so sets it to bool"""
        if self.label.lower() == "enable":
            self.value_type = "bool"
            self._check_toggle()
            
    def _check_toggle(self):
        """Checks to see if it should be checked depending on the value"""
        if self.value:
            self.checked = True
        else:
            self.checked = False
                
    def check_radio(self):
        """Checks the label to see if the value_type should be set to radio"""
        if self.label.lower() == "log_linear":
            self.value_type = "radio"
        if self.label.lower() == "linear_log":
            self.value_type = "radio"    
        
class ConfigModel(HasTraits):
    """This class contains the methods called by listView in the enaml file.

    This class the listModel variable which is sent to the listView in enanml
    and generates what is display in that list. The other variables in this
    class help control what data is displayed when an attribute is clicked on.
    """
    names,data,comments,labels = [],[],[],[]
    model = ListModel(names)
    index = 0
    current_name = Str
    current_data = List
    current_comment = Str
    current_label = Any
    fields = List
    key = 'display_defaults'
        
    def __init__(self):
        """Initializes 4 instances of FieldModel to be used by configform in enaml"""
        for i in range(4):
            self.fields.append(FieldModel())
        
    def set_file(self, json):
        """This method is called when a file is loaded and updates the listModel.
        
        This method takes the file selected and instantiates a copy of the
        json_import class. It uses that class to generate the list of names
        and data available and updates the listModel with that data.
        """
        self.loadedFile = json
        self.json_file = json_config(json)
        self.names = self.json_file.get_attrs()
        self.data = self.json_file.get_all_data()
        self.model = ListModel(self.names)
        
    def set_index(self, new):
        """This method is called when a new item is seleceted and updates values.
        
        When a new attribute is selected in the listModel, this method receives the
        new index and uses it to the update the current_name and current_data and uses
        the name to look up the labels for that attribute. This method also writes the
        data before it does all of this to avoid a race-condition.
        """
        self.write_data()
        self.index = new.row
        self.current_name = self.names[self.index]         
        self.current_data = self.data[self.index]
        self.labels = self.json_file.get_labels(self.current_name)
        self.labels = map(unicode.capitalize,self.labels)        
        self._populate_fields()       
    
    def _populate_fields(self):
        """This method populates all visible fields when a new attribute is selected
        
        This method is called by set_index, check to see how many fields need to be
        visible and updates the visible fields to the values of the new data. It
        the checks to see if these values are boolean types or require a radio
        button.
        """
        for i in range(len(self.current_data)):
            self.fields[i].set_value_type(self.type_check(self.current_data[i]))
            self.fields[i].set_visible(True)
            self.fields[i].set_value(self.current_data[i])
            self.fields[i].set_label(self.labels[i])
            self.fields[i].check_bool()
            self.fields[i].check_radio()
        for i in range(len(self.fields)-len(self.current_data)):
            self.fields[-(i+1)].set_visible(False)
            
    def type_check(self, value):
        """This methods returns the value_type of a given value.
        
        This method receives each value present for an attribute and returns
        their value type to be used by the validator.
        """
        value_type = type(value)
        if value_type == int:
            return "int"
        if value_type == float:
            return "float"
        if value_type == str:
            return "str"
        if value_type == unicode:
            return "unicode"
        
    def write_data(self):
        """This method updates the changes made to values to the list of data.
        
        This method is called right before changing to a new index and saves any
        changes made to the previous attribute to the data list.
        """
        
        for i in range(len(self.current_data)):
            self.data[self.index][i] = self.fields[i].value            
        
    def save_config(self):
        """This method saves the changes made to the configuration file.
        
        When the save button is pressed on the gui, this method is called,
        takes all the update values and uses them to update the dictionary
        and then saves the changes to the file.
        """
        self.write_data()
        currentPath,fileName = os.path.split(self.loadedFile)
        backupName = currentPath + '/' + os.path.splitext(fileName)[0] + '.bak'
        os.rename(self.loadedFile,backupName)
        attrs = self.json_file.get_attrs()
        for i in range(len(self.data)):
            labels = self.json_file.get_labels(attrs[i])
            for j in range(len(labels)):
                self.json_file.set_value(attrs[i],labels[j],self.data[i][j])
        data = json.dumps(self.json_file.obj, sort_keys=True, indent =4)
        f = open(self.loadedFile,'w')
        f.write(data)
        
    def save_config_as(self, path):
        """This method saves the changes made to the configuration file.
        
        When the save button is pressed on the gui, this method is called,
        takes all the update values and uses them to update the dictionary
        and then saves the changes to the file.
        """
        self.write_data()
        attrs = self.json_file.get_attrs()
        for i in range(len(self.data)):
            labels = self.json_file.get_labels(attrs[i])
            for j in range(len(labels)):
                self.json_file.set_value(attrs[i],labels[j],self.data[i][j])
        data = json.dumps(self.json_file.obj, sort_keys=True, indent =4)
        f = open(path,'w')
        f.write(data)
        f.close()
