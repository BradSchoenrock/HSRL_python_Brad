from datetime import datetime
from PyQt4.QtCore import *
from PyQt4.QtNetwork import *
from RawDataRay import RawDataRay
import sys
import time
import signal
import struct

#
# RawDataClient is a subclass of QThread which opens a socket connection
# to read from a server providing HSRL raw data (counts from the four HSRL 
# accumulators). When a new accumulation is available, RawDataClient emits a 
# newRawData(RawDataRay) signal.
#
# Simple usage will look something like:
#
#     def handlerMethod(ray):
#         assert isinstance(ray, RawDataRay), "ray is not a RawDataRay!"
#         print "Got a new ray!"
#
#     app = QCoreApplication(sys.argv)
#     readerThread = RawDataClient()
#     QObject.connect(readerThread, SIGNAL("newRawData(PyQt_PyObject)"), 
#                     handlerMethod)
#     readerThread.start()
#     sys.exit(app.exec_())
class RawDataClient(QThread):
    def __init__(self, host = "localhost", port = 44444):
        QThread.__init__(self)
        self.__rawServerHost = host
        self.__rawServerPort = port
        self.__socket = None
    
    def __connectToServer(self):
        if (self.__socket == None):
            self.__socket = QTcpSocket()
            
        # Try to connect until we're successful
        while (True):
            self.__socket.connectToHost(self.__rawServerHost, self.__rawServerPort);
            if (self.__socket.waitForConnected(1000)):
                # When connecting through a tunneled port, the QTcpSocket will
                # look connected even if there's no server active on the remote
                # side. A waitForReadyRead() on the socket will report 
                # RemoteHostClosedError when this is the case. Test for this,
                # and only return if we're *really* connected.
                reallyConnected = (self.__socket.waitForReadyRead(50) or 
                    self.__socket.error() <> QAbstractSocket.RemoteHostClosedError)
                if reallyConnected:
                    # We actually got connected!
                    return
            # Send a notification that we're waiting for a server to be available.
            print "waiting for raw data server to be available at", \
                self.__rawServerHost + ':' + str(self.__rawServerPort)
            self.sleep(5);
    
    def run(self):
        HEADER_LEN = 34 # size of the metadata in the version 1 raw data stream
        dataBuf = ""
        expectedLen = 0
        
        while (True):
            try:
                assert (self.__socket.state() == QAbstractSocket.ConnectedState)
            except:
                self.__connectToServer()
                
            # If no data are available, wait for up to five seconds for more to
            # show up, then try to reconnect.
            while (self.__socket.atEnd() and 
                   not self.__socket.waitForReadyRead(5000)):
                print "Reconnecting after no new data for a while"
                # Try to reconnect to the server
                self.__connectToServer()
            
            # Get what's available on the socket
            dataBuf += self.__socket.readAll()
            
            # Get more if we don't even have the header yet
            if (len(dataBuf) < HEADER_LEN):
                continue

            # If we haven't unpacked the header yet, do so now.
            if (expectedLen == 0):
                # Unpack the header of the raw data. All numeric values are
                # delivered in little-endian order.
                #    4 bytes - "HSRL"
                #    1 byte  - HSRL raw data version number
                #    4 bytes - unsigned number of gates
                #    8 bytes - unsigned ray time in microseconds
                #              since 1970-01-01 00:00:00 UTC
                #    4 bytes - unsigned number of good shots in the ray/accumulation
                #    4 bytes - unsigned transmitted energy monitor counts 
                #              (integrated over the ray/accumulation)
                #    4 bytes - floating point bin width, m
                #    4 bytes - floating point dwell time, s
                #    1 byte  - telescope pointing direction: 0->down, 1->up
                header = struct.unpack("<4sBLQLLffB", dataBuf[:HEADER_LEN])
                hsrl = header[0]
                version = header[1]
                nGates = header[2]
                rayTime = header[3]
                nGoodShots = header[4]
                xmitEnergyCounts = header[5]
                binWidth = header[6]
                dwellTime = header[7]
                pointingDir = header[8]
                
                # Validate that everything starts with "HSRL"
                if (hsrl == "HSRL"):
                    # Expected length is header size plus nGates of 4-byte data 
                    # for each of the four raw data channels (combined_hi, 
                    # combined_lo, cross, and molecular)
                    expectedLen = HEADER_LEN + (4 * nGates) * 4
                else:
                    print "'HSRL' leader not found, trying again..."
                    hsrlLoc = dataBuf.indexOf("HSRL")
                    if (hsrlLoc >= 0):
                        dataBuf = dataBuf[hsrlLoc:]
                    else:
                        dataBuf = ""
                    continue
                
                # We only read version 1 data!
                if (version != 1):
                    print 'Raw data are HSRL raw data version ' + version + \
                        ', not version 1 as expected!'
                    sys.exit(1)

                # Convert rayTime from microseconds since 1970-01-01 00:00:00 UTC
                # into a datetime
                rayTime = datetime.utcfromtimestamp(1.0e-6 * rayTime)

            # Go back for more if we don't have a full ray yet
            if (len(dataBuf) < expectedLen):
                continue
            
            # We now have the complete ray, so unpack the data arrays from it
            # The data arrays hold counts from the four accumulator channels of 
            # HSRL: "combined_hi", "combined_lo", "cross", and "molecular".
            # Each is an array of nGates 4-byte unsigned integer counts.
            format = "<" + str(nGates)  + "L"   # nGates 4-byte little-endian unsigned ints
            combined_hi = struct.unpack_from(format, dataBuf, HEADER_LEN);
            combined_lo = struct.unpack_from(format, dataBuf, HEADER_LEN + (4 * nGates))
            cross = struct.unpack_from(format, dataBuf, HEADER_LEN + 2 * (4 * nGates))
            molecular = struct.unpack_from(format, dataBuf, HEADER_LEN + 3 * (4 * nGates))
            
            
            # Construct a RawDataRay and emit the newRawDataRay(RawDataRay)
            # signal to pass it on to connected slots.
            newRay = RawDataRay(rayTime, nGates, nGoodShots, xmitEnergyCounts,
                                binWidth, dwellTime, pointingDir,
                                combined_hi, combined_lo, cross, molecular)
            self.emit(SIGNAL("newRawDataRay(PyQt_PyObject)"), newRay)
            
            # Keep leftover bytes if any, and go back for more
            dataBuf = dataBuf[expectedLen:]
            expectedLen = 0

if __name__ == "__main__":
    class Handler(QObject):
        def __init__(self, parent = None):
            QObject.__init__(self, parent)
        def printRayInfo(self, ray):
            assert isinstance(ray, RawDataRay), "ray is not a RawDataRay!"
            print str(ray.nBins()) + "-gate ray at " + str(ray.rayTime())
            
    # create a Qt core application
    app = QCoreApplication(sys.argv)
    
    # Allow ^C to shut down the application
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    # Create a thread to get the raw data
    rawDataClient = RawDataClient()
    
    # Create a handler object
    handler = Handler()

    # Connect incoming data to the handler's printRayInfo slot
    QObject.connect(rawDataClient, SIGNAL("newRawDataRay(PyQt_PyObject)"),
                    handler.printRayInfo)
    
    # Start the raw data client thread
    rawDataClient.start()
    
    # start the event loop
    sys.exit(app.exec_())
