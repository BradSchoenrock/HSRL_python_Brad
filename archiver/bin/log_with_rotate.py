#!/usr/bin/python
import sys
import logging 
import logging.handlers

def log_with_rotate(log_filename, maxBytes=200*1000*1000):
    # set up a specific logger with our desired output level
    my_logger = logging.getLogger('MyLogger')
    my_logger.setLevel(logging.DEBUG)

    # Add the log message handler to the logger
    handler = logging.handlers.RotatingFileHandler(
              log_filename, maxBytes=maxBytes, backupCount=5)

    my_logger.addHandler(handler)
    while 1:
	line = sys.stdin.readline()
	if len(line) == 0: break
	my_logger.debug(line[0:-1])   # remove trailing newline


if __name__ == '__main__':
    if len(sys.argv) != 2:
	print 'usage: %s log_filename ' % sys.argv[0]
	sys.exit(1)
    log_with_rotate(sys.argv[1])


