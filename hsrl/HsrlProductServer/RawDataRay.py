#!/usr/bin/python
from datetime import datetime
from copy import copy
from PyQt4.QtCore import *

class RawDataRay(QObject):
    def __init__(self, rayTime, nBins, nGoodShots, xmittedEnergyCounts,
                 binWidth, dwellTime, pointingDir, 
                 combinedHiData, combinedLoData, crossData, molecularData):
        """
        Construct a RawDataRay from the given arguments:
        
        rayTime                - datetime object containing the ray time
        nBins                  - number of bins in the data
        nGoodShots             - number of good laser shots in the accumulation
        xmittedEnergyCounts    - raw counts of transmitted energy summed
                                 for all shots in the ray
        binWidth               - width of a data bin, m
        dwellTime              - period over which the accumulation was sampled
        pointingDir            - telescope pointing direction (0->down, 1->up)
        combinedHiData         - nBins of raw counts from the 'combined_hi' accumulator channel
        combinedLoData         - nBins of raw counts from the 'combined_lo' accumulator channel
        crossData              - nBins of raw counts from the 'cross' accumulator channel
        molecularData          - nBins of raw counts from the 'molecular' accumulator channel
        """

        QObject.__init__(self)
        if (not isinstance(rayTime, datetime)):
            raise TypeError('rayTime must be of type datetime')
        self.__rayTime = copy(rayTime)
        self.__nBins = nBins
        self.__nGoodShots = nGoodShots
        self.__xmittedEnergyCounts = xmittedEnergyCounts
        self.__binWidth = binWidth
        self.__dwellTime = dwellTime
        self.__pointingDir = pointingDir
        # Copy the data arrays, validating that the length of each array matches
        # the given gate count
        assert (len(combinedHiData) == nBins), "bad length for combined_hi data"
        self.__combinedHiData = combinedHiData
        assert (len(combinedLoData) == nBins), "bad length for combined_lo data"
        self.__combinedLoData = combinedLoData
        assert (len(crossData) == nBins), "bad length for cross data"
        self.__crossData = crossData
        assert (len(molecularData) == nBins), "bad length for molecular data"
        self.__molecularData = molecularData
    
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
    
    def nGoodShots(self):
        """
        The number of good laser shots accumulated in this ray
        """
        return self.__nGoodShots
    
    def xmittedEnergyCounts(self):
        """
        The raw counts from the transmitted energy monitor summed over all 
        laser shots included in this ray. 
        """
        return self.__xmittedEnergyCounts
    
    def binWidth(self):
        """
        The bin width for this ray, m
        """
        return self.__binWidth
    
    def dwellTime(self):
        """
        The dwell time for this ray, s
        """
        return self.__dwellTime
    
    def pointingDir(self):
        """
        Pointing direction for this ray (0->down, 1->up)
        """
        return self.__pointingDir
    
    def combinedHiData(self):
        """
        Return the list of counts from the 'combined_hi' channel for this
        ray.
        """
        return self.__combinedHiData
    
    def combinedLoData(self):
        """
        Return the list of counts from the 'combined_lo' channel for this
        ray.
        """
        return self.__combinedLoData
    
    def crossData(self):
        """
        Return the list of counts from the 'cross' channel for this
        ray.
        """
        return self.__crossData
    
    def molecularData(self):
        """
        Return the list of counts from the 'molecular' channel for this
        ray.
        """
        return self.__molecularData
