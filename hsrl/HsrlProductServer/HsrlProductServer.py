#!/usr/bin/env python
import sys
from optparse import OptionParser
from PyQt4.QtCore import *
from RawDataClient import RawDataClient
from HsrlProductGenerator import HsrlProductGenerator
from MultiClientTcpServer import MultiClientTcpServer
import struct
from time import mktime
import signal

class StopAppOnCount(QObject):
    def __init__(self, app, stopCount):
        QObject.__init__(self)
        self.__app = app
        self.__stopCount = stopCount
        self.__count = 0
        
    def incrementCount(self):
        self.__count = self.__count + 1
        if (self.__count == self.__stopCount):
            self.__app.exit(0)

def publishProductRay(pRay):
    """
    Most of the real work happens here. The incoming ProductRay is converted 
    into on-the-wire form described below. Then the on-the-wire form is 
    published via our server.

    
    Over-the-wire product format
    ============================
    All numeric values in the on-the-wire format are delivered in little-endian 
    order.
    
    8 bytes        - char: 'HSRLPROD'
    8 bytes        - unsigned int: time in microseconds since 1970-01-01 
                     00:00:00 UTC
    4 bytes        - unsigned int: number of bins in data
    4 bytes        - float: bin width, m
    4 bytes        - float: dwell time, s
    4 bytes        - unsigned int: nProducts
    
    followed by nProducts product segments containing:
    
    32 bytes        - char: null-terminated string of product name
    32 bytes        - char: null-terminated string of product units
    nBins * 4 bytes - float[nBins]: data for the product
    """
    # Convert ray's datetime into microseconds since 1970-01-01 00:00 UTC
    secs = mktime(pRay.rayTime().timetuple()) + 1.0e-6 * pRay.rayTime().microsecond
    usecs = long(1.0e6 * secs) # 8-byte integer microseconds
    
    # Construct the header of the outgoing packet
    packet = ""
    packet += struct.pack("<8sQLffL", "HSRLPROD", usecs, pRay.nBins(), 
                          pRay.binSpacing(), pRay.dwellTime(),
                          pRay.nProducts())
    # Append the products
    for productname in pRay.productNames():
        product = pRay.product(productname)
        packet += struct.pack("31sc31sc", product.name(), '\0', 
                              product.units(), '\0')

        format = "<" + str(pRay.nBins())  + "f"   # nBins little-endian 4-byte floats
        packet += struct.pack(format, *product.data())
        
    # Finally, send out the packet via our server
    server.sendToClients(packet)

# Parse our command line options
parser = OptionParser()

defRawSource = "localhost:44444"
defOutputPort = 55555
defNTimeAvg = 1
defNBinAvg = 1
defMaxBin = 4000
parser.set_defaults(rawDataSource = defRawSource, 
                    outputPort = defOutputPort,
                    nTimeAvg = defNTimeAvg,
                    nBinAvg = defNBinAvg,
                    maxBin = defMaxBin)

parser.add_option("--rawDataSource", type = "string", dest = "rawDataSource", 
                  help = "server:port providing raw data, default " + defRawSource)

parser.add_option("--outputPort", type = "int", dest = "outputPort",  
                  help = "output port for products, default " + str(defOutputPort))

parser.add_option("--nTimeAvg", type = "int", dest = "nTimeAvg", 
                  help = "number of raw rays per time average, default " + str(defNTimeAvg))

parser.add_option("--nBinAvg", type = "int", dest = "nBinAvg", 
                  help = "number of raw bins per range average, default " + str(defNBinAvg))

parser.add_option("--maxBin", type = "int", dest = "maxBin",
                  help = "maximum number of raw bins to process, default " + str(defMaxBin))


(options, args) = parser.parse_args()

# Parse the raw data host and port number from options.rawDataSource
(rawHost, rawPort) = tuple(options.rawDataSource.split(":"))
try:
    rawPort = int(rawPort)
except ValueError:
    print "Bad raw data port '" + rawPort + "'"
    sys.exit(1)

# Get output port, nTimeAvg, and nBinAvg from options    
outputPort = options.outputPort
nTimeAvg = options.nTimeAvg
nBinAvg = options.nBinAvg
maxBin = options.maxBin

# create a Qt core application, passing it the remaining command line arguments
app = QCoreApplication(args)

# Client to read raw data
reader = RawDataClient(rawHost, rawPort)

# Product generator. 
prodGen = HsrlProductGenerator(maxBin, nTimeAvg, nBinAvg)

# create the server, which will send the data to clients on the net
server = MultiClientTcpServer(outputPort)

# Send raw data to the product generator. Use a QueuedConnection to make
# the processing happen in the prodGen thread!
QObject.connect(reader, SIGNAL("newRawDataRay(PyQt_PyObject)"), 
                prodGen.handleRawRay, Qt.QueuedConnection)

# Pass products from prodGen to publishProductRay()
QObject.connect(prodGen, SIGNAL("newProductRay(PyQt_PyObject)"), 
                publishProductRay, Qt.QueuedConnection)

## Stop app after a set number of product rays are generated
#appStopper = StopAppOnCount(app, 25)
#QObject.connect(prodGen, SIGNAL("newProductRay(PyQt_PyObject)"),
#                appStopper.incrementCount, Qt.QueuedConnection)

# Allow ^C to shut down the application
signal.signal(signal.SIGINT, signal.SIG_DFL)

# start the raw data reader and product generator
reader.start()

# start the event loop
sys.exit(app.exec_())
