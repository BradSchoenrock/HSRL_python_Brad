#!/usr/bin/python
from datetime import datetime
from copy import copy
from PyQt4.QtCore import *


class Product():
    """
    Product is a simple class supporting a data "product" with name, units,
    and data.
    """
    def __init__(self, name, units, data):
        """
        Construct a product using the given name string, units string, and
        data. The data provided must support len() and [] operations.
        """
        self.__name = name
        self.__units = units
        # validate that the len() and [] operations work with the given 
        # data
        try:
            nbins = len(data)
            val = data[-1]
        except:
            raise TypeError("The data provided for Product must support len() " +
                            "and [] operations.")
        self.__data = data
        
    def name(self):
        """
        Return the name of the product
        """
        return self.__name
    
    def units(self):
        """
        Return the units of the product
        """
        return self.__units
    
    def data(self):
        """
        Return the data for the product, which will support indexed access
        """
        return self.__data
    
    def nBins(self):
        """
        Return the number of bins in the data for this product.
        """
        return len(self.__data)
        
    

class ProductRay(QObject):
    def __init__(self, rayTime, nBins, binSpacing, dwellTime):
        """
        Construct a ProductRay from the given ray time (in datetime 
        form), gate/bin count, bin spacing in meters, and dwell time in 
        seconds. The ray starts out with no associated products, but they can 
        be added using the addProduct() method.
        """
        QObject.__init__(self)
        if (not isinstance(rayTime, datetime)):
            raise TypeError('rayTime must be of type datetime')
        self.__rayTime = copy(rayTime)
        self.__nBins = nBins
        self.__binSpacing = binSpacing
        self.__dwellTime = dwellTime
        self.__products = {}

    def rayTime(self):
        """
        Return the datetime for this ray
        """
        return self.__rayTime
    
    def nBins(self):
        """
        Return the number of bins in this ray
        """
        return self.__nBins
    
    def binSpacing(self):
        """
        Return the bin spacing for the ray, m
        """
        return self.__binSpacing
        
    def dwellTime(self):
        """
        Return the dwell time of the ray, s
        """
        return self.__dwellTime
        
    def addProduct(self, product):
        # Verify that product is of type Product, and that its bin count 
        # matches ours.
        assert isinstance(product, Product), "'product' must be of type Product"
        assert (product.nBins() == self.__nBins), \
            "In ProductRay.addProduct(): product.nBins() must match ProductRay.nBins()"
        # Add the product to our dictionary of products
        self.__products[product.name()] = product
        
    def productNames(self):
        """
        Return a list containing the names of all our products
        """
        return sorted(self.__products.keys())
    
    def nProducts(self):
        """
        Return the number of products in this ray
        """
        return(len(self.__products.keys()))
    
    def product(self, name):
        """
        Return a Product object for the named product.
        """
        return self.__products[name]
