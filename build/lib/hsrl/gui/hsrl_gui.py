# HSRL GUI 
import os
from traits.api import HasTraits, Range, on_trait_change, Str, Date, Int
from hsrl.data_stream.rti import Rti
from lg_base.core.getTime import get_time
from lg_base.core.parsers import parse_dates, get_workingdirs
from lg_base.core.configModel import ConfigModel

# WARNING - because of the class level 'get_workingdirs()' call, 
# HSRL_CONFIG and HSRL_INSTRUMENT
# must be set before calling this routine.
#
# TODO - revisit this, and set the class static variables later
    
class HSRL_Model(HasTraits):
    """Model class which is used by enaml to snyc the GUI to the commandline.
    
    This class contains the variables and functions which are called and accessed
    by the enaml frontend of the GUI and then sent to Rti after being called.
    Each variable represents a value that can be set/changed on the GUI and is 
    automatically updated whenever a value is changed on the GUI.
    """
    #These assignments needed by enaml.
    working_dir = get_workingdirs()
    lidar = Str
    plot_length = Range(low=0.0)
    min_altitude = Range(low=-0.5)
    max_altitude = Range(low=-0.5)
    optical_depth = Range(low=-0.5)
    dates = parse_dates()
    year_array = dates[0]
    year = year_array[0]
    month_array = dates[1]
    month = month_array[0][0]
    day_array = dates[2]
    day = day_array[0][0][0]
    date_time = get_time(year, month, day, working_dir)
    filename = Str
    
    def update_time(self):
        """Updates the earliest available time.
        
        This module is called after year, month or day are changed and 
        calls get_time to retrieve the earliest avaiable time of data
        for the new date.
        """
        self.date_time = get_time(self.year, self.month, 
                                  self.day, self.working_dir)
        return self.date_time
    
    def start_plot(self):
        """Starts the plotting routine when called
        
        Sends the selected data from the GUI to an instance of Rti
        and begins generating the plots when called
        """
        
        print 'start plot called:'
        #self.r=Rti('gvhsrl','18-Aug-11 00:00', 2, 0, 15)
        self.r=Rti(self.lidar, self.date_time.strftime("%d-%b-%y %H:%M"), 
              self.plot_length, self.min_altitude, self.max_altitude, 
              self.optical_depth, display=self.filename)
    
    def next_plot(self):
        """Updates the generated plots to the next time inverval
        
        When called after plots are generated, it advances the 
        current plots forward equal to the time interval
        """
        
        self.r.next()
             
    
def start_gui():
    """Generates the HSRL_GUI when called.
    
    On call, generates the HSRL_GUI and check to see if the environment
    variables HSRL_CONFIG and HSRL_INSTRUMENT are set.
    Afterwards, it generates an instance of HSRL_MODEL which contains
    the variables which will be passed to Rti. It will then generate
    an instance of ConfigModel which controlls config file editing.
    """
    
    if not os.getenv('HSRL_CONFIG') and not os.getenv('HSRL_INSTRUMENT'):
        print "HSRL_CONFIG and HSRL_INSTRUMENT environment variables are not set."
        return
    if not os.getenv('HSRL_CONFIG'):
        print "HSRL_CONFIG environment variable is not set."
        return
    if not os.getenv('HSRL_INSTRUMENT'):
        print "HSRL_INSTRUMENT environment variable is not set."
        return
        
    model = HSRL_Model(lidar='gvhsrl',plot_length=1,min_altitude=5,
                   max_altitude=20,optical_depth=10,
                   date_time=HSRL_Model.date_time)
    import enaml
    with enaml.imports():
        from hsrl_gui_view import HSRL_GUI

    msg = 'HSRL Plot Control'
    list_model = ConfigModel()
    window = HSRL_GUI(model=model,message=msg,list_model=list_model)
    window.show()

if __name__ == '__main__':   
    start_gui()
