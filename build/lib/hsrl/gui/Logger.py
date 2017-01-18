from traits.api import HasTraits, Str, Range, on_trait_change

#from enaml.stdlib.sessions import simple_app
#from enaml.qt.qt_local_server import QtLocalServer


class Logger(HasTraits):
    """
        class to log messages to an enaml TextEdit widget
    """    
    theMsg = Str    
    
    @on_trait_change('theMsg')
    def debug_print(self):
        """
        prints out a debug message when 'theMsg' changes        
        """
        print 'theMsg =', self.theMsg
    
    def doLog(self, message):
        self.theMsg = message        

if __name__ == '__main__':
    import enaml    
    with enaml.imports():
        from log_view import LogView
    
    logger = Logger(theMsg = 'test message')    
    
    window = LogView(logger = logger)
    window.show()
    
