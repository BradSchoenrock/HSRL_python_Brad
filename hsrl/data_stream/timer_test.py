from threading import Timer
import os

# create a pipe which we load with NLs 
rfd=None
wfd=None
rpipe=None
wpipe=None


def keybanger():
   global timmy   
   wpipe.write('\n')
   wpipe.flush()
   print ':'
   timmy = Timer(5.0, keybanger)   # may not need this 
   timmy.start()

if __name__=='__main__':
	rfd, wfd = os.pipe()
	rpipe, wpipe = os.fdopen(rfd, 'r'), os.fdopen( wfd, 'w' )
	timmy = Timer(5.0, keybanger)
	timmy.start()

# then replace raw_input with 
 #while True: 
    # ...
 #   cmd = rpipe.read(1) 
 #   print '.'
