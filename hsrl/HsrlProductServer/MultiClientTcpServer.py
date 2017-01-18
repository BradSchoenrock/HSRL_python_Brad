import sys
import signal
from PyQt4.QtCore import *
from PyQt4.QtNetwork import *

class MultiClientTcpServer(QTcpServer):
    """
    MultiClientTcpServer delivers TCP packet data to all clients which
    connect to its server port. The port number for data delivery is 
    specified when the server is created.
    """
    def __init__(self, port, parent = None):
        QTcpServer.__init__(self, parent)
        # Listen for connections on the given port
        if (not self.listen(QHostAddress(QHostAddress.Any), port)):
            print "Failure to listen on port " + str(port) + ":", self.errorString()
            sys.exit(1)
        else:
            print "Listening on port " + str(port)

        # Start with an empty list of clients
        self.__clients = []
        
        # Add clients to our list when a newConnection() signal is emitted
        QObject.connect(self, SIGNAL("newConnection()"), self.__addNewClients)
        
    @classmethod
    def __clientName(cls, client):
        name = client.peerAddress().toString() + ":" + str(client.peerPort())
        return name
        
    def __addNewClients(self):
        """
        Slot to add all pending clients to our client list
        """
        while True:
            client = self.nextPendingConnection()
            if (client == None):
                break
            
            # Add this socket to our list of clients
            self.__clients.append(client);
            
            # When the client disconnects, remove it from our list of clients.
            QObject.connect(client, SIGNAL("disconnected()"), self.__removeClient)

            print "connection from", self.__clientName(client)
            
    def __removeClient(self):
        """
        Slot to remove the signal sender from our client list.
        """
        client = self.sender()
        if (client in self.__clients):
            self.__clients.remove(client)
            
            print "disconnect from", self.__clientName(client)
            
    def sendToClients(self, data):
        """
        Send the given data, of type QByteArray, to each of the connected 
        clients.
        """
        for client in self.__clients:
            result = client.write(data)
            if (result < 0):
                print "Error writing to", self.__clientName(client), "-", client.errorString()
            elif (result <> len(data)):
                print "Only wrote", result, "of", len(data), "bytes to", self.__clientName(client)
        
if __name__ == "__main__":
    # create a Qt core application
    app = QCoreApplication(sys.argv)
    
    # Allow ^C to shut down the application
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    # create the server
    server = MultiClientTcpServer(44444)
    
    # start the event loop
    sys.exit(app.exec_())
    
