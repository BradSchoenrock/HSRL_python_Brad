import sys
import traceback
sys.path.append('../data_stream')
sys.path.append('../calibration')

import thread
import numpy
from copy import copy
from PyQt4.QtCore import *
from PyQt4.QtNetwork import *
from datetime import datetime
from matplotlib.dates import date2num
import hsrl.data_stream.hsrl_read_utilities as hru
import hsrl.calibration.calibration_utilities as cu
import lg_base.formats.calvals as cru
import sounding_utilities as su
import processing_utilities as pu
import tz_utilities as tzu
from ProductRay import Product
from ProductRay import ProductRay

SPEED_OF_LIGHT = 299792458.     # m/s

# cal_vectors is ripped from main_routines.py
class cal_vectors(object):
    def __init__(self, instrument, time, max_range_bin):
        """initialize calibration vector structure"""
        self.instrument = instrument
        self.max_range_bin = max_range_bin
        self.read(time)          


    def read(self, time):
        """read baseline, geo, diff_geo, and i2_scan files"""
        self.baseline = hru.read_baseline(self.instrument, time)
        self.baseline.data = self.baseline.data[:self.max_range_bin, :]
        self.geo = hru.read_geo_corr(self.instrument, time)
        if self.geo.data!=None:
            self.geo.data = self.geo.data[:self.max_range_bin, :]
        self.diff_geo = hru.read_diff_geo(self.instrument, time)
        if self.diff_geo.data!=None:
            self.diff_geo.data = self.diff_geo.data[:self.max_range_bin, :]
        self.i2scan = hru.read_i2_scan(self.instrument, time)
        # first expiration time
        self.expire_time = min(self.baseline.expire_time,
                               self.geo.expire_time,
                               self.diff_geo.expire_time,
                               self.i2scan.expire_time)

class RawRayMatchException(Exception):
    '''
    RawRayMatchException can be raised if there is a mismatch in raw rays
    which does not allow them to be processed together.
    '''
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ProductGeneratorThread(QThread):
    def __init__(self, instrument, nBins, nTimeAvg, nRangeAvg):
        '''
        Instantiate a HsrlProductGeneratorThread.

        instrument - the name of the HSRL instrument providing raw data
                     (e.g., 'gvhsrl')
        nBins - the number of bins to process from the raw data
        nTimeAvg - the number of raw rays to average in time when processing
        nRangeAvg - the number of bins to average in range when processing
        '''
        QThread.__init__(self)

        # Mutex for data protection between threads, and condition to allow
        # for waking processing
        self.__mutex = QMutex()
        self.__dataReady = QWaitCondition()

        # Initialize 
        self.__rawRays = None
        self.__nBins = nBins
        self.__nTimeAvg = nTimeAvg
        self.__nRangeAvg = nRangeAvg
        self.__instrument = copy(instrument)

        # We cannot load calibration information until we have a "real" data
        # time to work with. Just create placeholders for now.
        self.__calvals = None
        self.__calConstants = None
        self.__currentCal = None
        self.__soundings = None
        self.__currentSounding = None
        self.__calibration = None
        self.__inversionCoefs = None

        # Processing selection array. This will be created in 
        # __initCalibrations()
        self.__psel = None

        # KLUGE-and-a-half
        # Even though we aren't dealing with displays here, we still have
        # "display defaults" in order to use processing_utilities.process_data().
        # Try to build a minimal subset here that will make process_data() 
        # happy.
        self.__klugeDisplayDefaults = {}
        self.__klugeDisplayDefaults['mask_image'] = [0, 250]
        self.__klugeDisplayDefaults['show_cal_pulse'] = None

        # Mutex for data protection between threads, and condition to allow
        # for waking processing
        self.__mutex = QMutex()
        self.__startProcessing = QWaitCondition()

    def __del__(self):
        # Hmmm, this destructor doesn't really work...
        if (self.__workThread):
            self.__workThread.exit(0)

    def generateProductRay(self, rawRays):
        '''
        Slot which accepts an incoming list of RawDataRay-s of raw accumulator 
        data and processes them in the working thread. A newProductRay(ProductRay)
        signal is emitted when the thread completes the processing.
        '''
        QMutexLocker(self.__mutex)
        self.__rawRays = rawRays
        self.__dataReady.wakeAll()

    def isBusy(self):
        '''
        Returns true iff this thread is busy generating a product ray.
        '''
        # If we can get a lock on the mutex, we're not currently busy
        if (self.__mutex.tryLock()):
            self.__mutex.unlock()
            return False
        else:
            return True

    def __loadCalForTime(self, time):
        """
        Load calibration/sounding information valid for the given time,
        in the matplotlib form (1 + days since 0001-01-01 00:00:00 UTC).
        """
        # Initialize if necessary, using the given time as the start time
        if (self.__calvals == None):
            self.__initCalibrations(time)
        # Create a bogus raw data time_z_data array, which is required but will 
        # go unused in our call to update_calibrations()...
        bogusRawData = tzu.time_z_data()
        # Update to get the calibration for the given time
        calDict = cu.update_calibrations(self.__instrument, time, 1,
                                         bogusRawData,
                                         self.__calConstants, self.__soundings,
                                         self.__currentCal, self.__calvals,
                                         self.__currentSounding)
        # Pull the updated bits from the return dictionary
        self.__calConstants = calDict['rs_constants']
        self.__currentCal = calDict['rs_cal']
        self.__currentSounding = calDict['sounding']

    def __initCalibrations(self, time):
        """
        Initialize our calibration information with the given time as our 
        starting time. Time should be in matplotlib form (1 + days since 
        0001-01-01 00:00:00 UTC).
        """
        # read calvals file of system calibration constants
        self.__calvals = cru.cal_file_reader(self.__instrument)

        # Initialize system constants at the given time
        self.__calConstants = \
            cru.select_sys_constants(self.__calvals, time)

        # KLUGE: select_sys_constants() does not provide extinction_bin_delta 
        # (which comes from display_defaults.json/all_plots.json), but the 
        # value is required by processing_utilities.compute_extinction(). We 
        # don't load either of those .json files, so just put in a fixed value 
        # for now.
        self.__calConstants['extinction_bin_delta'] = 2

        # The binwidth in the constants is in seconds. Convert to meters.
        binwidth = 0.5 * SPEED_OF_LIGHT * self.__calConstants['binwidth']

        # Calculate our range resolution based on our range averaging factor
        self.__rangeRes = self.__nRangeAvg * binwidth

        # Read the soundings file and extend the soundings to 50 km if needed 
        # using climatology.
        # Soundings range from one before the start time to the end of file.
        # Values returned with altitude resolution set to bin spacing.
        args = []
        args.append(self.__instrument)
        args.append(self.__calConstants['sounding_type'])
        args.append(self.__calConstants['default_raob_id'])
        args.append(time)
        # It is IMPERATIVE that the altitude resolution given here matches
        # the one we will pass to processing_utilities.process_data(), since
        # the dimensions of the soundings built here are assumed to match what
        # will be generated there.
        args.append(self.__rangeRes)
        args.append(50000)   # extend soundings up to 50 km
        try:
            self.__soundings = su.read_sounding_file(args)
        except:
            print "Error reading sounding file:"
            traceback.print_exc()
            sys.exit(1)

        # Get the specific sounding for the given time
        self.__currentSounding = \
            tzu.select_time_z_data_object(self.__soundings, time)

        # Initialize our calibration
        self.__currentCal = cal_vectors(self.__instrument, time, self.__nBins)

        # Processing selection array. Start with all enabled and and adjust 
        # below.
        self.__psel = numpy.ones(9)

        # use the ratio of combined to molecular channel sensitivity from our cal
        self.__psel[4] = self.__calConstants['i2_scan_adjustment']
        # turn signal_in_dark_count correction off; it's not currently working
        self.__psel[6] = 0

        # Use witschas spectral model
        spectral_model = 'witschas'

        # Inversion coefficients are generated at the requested altitude 
        # resolution up to 50km
        self.__inversionCoefs = \
            cu.quick_cal(self.__currentCal.i2scan.data,
                         self.__currentCal.i2scan.constant, 
                         self.__currentSounding, 
                         self.__calConstants['wavelength'], 
                         spectral_model,
                         self.__psel)

    def __prf(self):
        """
        Return the lidar pulse repetition frequency, Hz. 
        This function returns 0.0 if no calibration has been loaded yet.
        """
        if (self.__calConstants is None):
            return 0.0
        else:
            return self.__calConstants['laser_rep_rate']

    def run(self):
        """
        Loop forever waiting for the __dataReady wait condition, then calling
        __processRawRays().
        """
        self.__mutex.lock()
        while True:
            # Wait until new data are ready, then process the raw rays
            self.__dataReady.wait(self.__mutex)
            try:
                # Generate a new ProductRay from the raw data rays we were 
                # given, and emit the newProductRay(ProductRay) signal
                pRay = self.__generateProductRay()
                self.emit(SIGNAL("newProductRay(PyQt_PyObject)"), pRay)
            except RawRayMatchException as rrme:
                print "Product ray not processed at", \
                      self.__rawRays[-1].rayTime(), "because: ", rrme
            except:
                print ""
                print "__generateProductRay() problem at", self.__rawRays[-1].rayTime()
                traceback.print_exc()
                sys.exit(1)

    def __generateProductRay(self):
        '''
        Process all raw rays in the self.__rawRays list and return a single 
        ProductRay.
        If all raw rays do not point the same direction, a RawRayMatchException
        is raised.
        '''
        rawRays = self.__rawRays
        nTimes = len(rawRays)
        nBins = self.__nBins

        # We need a time_z_data struct for the raw data to be sent to 
        # process_data(). Start with an empty one.
        rawData = tzu.time_z_data()

        # Add the nTimes x 1 arrays
        timeArrayShape = (nTimes, 1)
        rawData.times = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.GPS_MSL_Alt = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.telescope_pointing = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.seeded_shots = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.transmitted_energy = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.roll_angle = tzu.t_array.from_array(numpy.zeros(timeArrayShape))
        rawData.pitch_angle = tzu.t_array.from_array(numpy.zeros(timeArrayShape))

        # The raw count data go in nTimes x nBins arrays
        dataArrayShape = (nTimes, nBins)
        rawData.combined_hi_counts = tzu.tz_array.from_array(numpy.zeros(dataArrayShape))
        rawData.combined_lo_counts = tzu.tz_array.from_array(numpy.zeros(dataArrayShape))
        rawData.cross_pol_counts = tzu.tz_array.from_array(numpy.zeros(dataArrayShape))
        rawData.molecular_counts = tzu.tz_array.from_array(numpy.zeros(dataArrayShape))

        # Create an array of bin distances from the lidar (just a list of
        # bin number * bin width)
        binRanges = [b * rawRays[0].binWidth() for b in range(nBins)]

        # Add the bin altitudes array, which must be dimensioned (1 x nBins)
        # We enclose the 1-d list in [] to make it 2-d, and use that to 
        # instantiate a numpy.array().
        rawData.altitudes = tzu.z_array.from_array(numpy.array([binRanges]))

        # Set the processing start and end times to the first and last times
        # in the raw data. These times are kept in matplotlib format (number 
        # of days since 0001-01-01 00:00:00 UTC), since that's what is expected
        # by process_data()
        startTime = datetime(2500, 1, 1)    # 2500-01-01, i.e., way in the future
        endTime = datetime(1900, 1, 1)      # 1900-01-01, i.e., way in the past

        for i in range(nTimes):
            startTime = min(startTime, rawRays[i].rayTime())
            endTime = max(endTime, rawRays[i].rayTime())

        mplStartTime = date2num(startTime)
        mplEndTime = date2num(endTime)   

        # Load the calibration for the end of the processing period
        self.__loadCalForTime(mplEndTime)

        # Put necessary pieces from each of the raw rays into the rawData
        # struct which will be passed to process_data()
        for i in range(nTimes):
            # Add the time for this raw ray. Times must be represented in 
            # matplotlib (mpl) form: (1 + days since 0001-01-01 00:00:00 UTC), 
            # as generated by matplotlib.dates.date2num().
            mplTime = date2num(rawRays[i].rayTime())
            rawData.times[i, 0] = mplTime

            # Lidar altitude
            # @TODO: use actual aircraft altitude here
            rawData.GPS_MSL_Alt[i, 0] = self.__calConstants['lidar_altitude']

            # Make sure all raw rays are pointing the same way, or raise an
            # exception.
            if (rawRays[i].pointingDir() != rawRays[0].pointingDir()):
                raise RawRayMatchException('Not all raw rays to point in the same direction')
            # Telescope direction is dimensioned (nTimes x 1)
            # 0->down, 1->up
            rawData.telescope_pointing[i, 0] = rawRays[i].pointingDir()

            # Number of good laser shots accumulated in the raw ray
            rawData.seeded_shots[i, 0] = rawRays[i].nGoodShots()

            # Total transmitted energy from all laser shots accumulated in the
            # raw ray (in raw counts from the energy monitor)
            rawData.transmitted_energy[i, 0] = rawRays[i].xmittedEnergyCounts()

            # Aircraft pitch and roll angles, also dimensioned (nTimes x 1)
            # @TODO put real values here
            rawData.roll_angle[i, 0] = 0.0
            rawData.pitch_angle[i, 0] = 0.0

            # Add count data from the four accumulators as (1 x nBins) numpy 
            # arrays.
            rawData.combined_hi_counts[i] = numpy.array(rawRays[i].combinedHiData()[:nBins])
            rawData.combined_lo_counts[i] = numpy.array(rawRays[i].combinedLoData()[:nBins])
            rawData.cross_pol_counts[i] = numpy.array(rawRays[i].crossData()[:nBins])
            rawData.molecular_counts[i] = numpy.array(rawRays[i].molecularData()[:nBins])

        # Compute the time resolution and range resolution we need to pass
        # to process_data()
        if nTimes == 1:
            timeRes = 0.0   # process_data() magic value for no time averaging
        else:
            timeRes = rawRays[-1].dwellTime() * nTimes    # seconds

        # Get min and max altitudes of the data
        maxRange = nBins * rawRays[0].binWidth()
        lidar_alt = rawData.GPS_MSL_Alt[0, 0]
        if (rawData.telescope_pointing[0, 0] == 0):
            # Downward pointing
            min_alt = lidar_alt - maxRange
            max_alt = lidar_alt
        else:
            # Upward pointing
            min_alt = lidar_alt
            max_alt = lidar_alt + maxRange

        # Process the raw data
        with numpy.errstate(divide='ignore'):
            processedData = \
                pu.process_data(self.__instrument, mplStartTime, mplEndTime, 
                                min_alt, max_alt, 
                                timeRes, self.__rangeRes, 
                                rawData, self.__inversionCoefs, 
                                self.__calConstants, self.__currentCal, 
                                self.__currentSounding,
                                self.__psel, [], False, 
                                self.__klugeDisplayDefaults)

        # If we had to use zero time resolution above, get the *real* time
        # resolution now
        if timeRes == 0.0:
            timeRes = rawRays[-1].dwellTime()

        # Create a new ProductRay
        nProcessedBins = processedData.rs_inv.beta_a_backscat_par.shape[1]
        pRay = ProductRay(rawRays[-1].rayTime(), nProcessedBins, 
                          self.__rangeRes, timeRes)

        # Add aerosol backscatter cross-section to the ProductRay
        abcsData = processedData.rs_inv.beta_a_backscat_par[0] * \
            (1 + processedData.rs_inv.circular_depol[0])
        product = Product("aerosol backscatter cross-section", "m-1 sr-1",
                          abcsData)
        pRay.addProduct(product)

        # Add linear depolarization ratio to the ProductRay
        # From ../data_stream/display_utilities.py:
        #     linearDepol = circularDepol / (2 + circularDepol)
        circDepolData = processedData.rs_inv.circular_depol[0]
        linDepolData = circDepolData / (2 + circDepolData)
        product = Product("linear depolarization", "", linDepolData)
        pRay.addProduct(product)

        # Add the combined_hi accumulator channel data. KLUGE: Note that we 
        # will *not* generally have a bin-to-bin correspondence with the 
        # processed data, since the processed data exists from min_alt to
        # max_alt, while the raw data remains in raw data space.
        combinedHiData = [float(i) for i in rawData.combined_hi_counts[0,0:nProcessedBins]]
        product = Product("combined_hi (raw)", "counts", combinedHiData)
        pRay.addProduct(product)

        return pRay

class HsrlProductGenerator(QObject):
    """
    HsrlProductGenerator accepts incoming raw HSRL accumulator data in the form 
    of RawDataRay objects delivered to its processRawRay() socket. Each time it 
    accumulates nTimeAvg raw rays, it makes a new ProductRay and emits it via 
    its newProductRay(ProductRay) signal.
    """

    def __init__(self, maxBins, nTimeAvg, nRangeAvg, instrument = 'gvhsrl'):
        """
        Instantiate a HsrlProductGenerator.

        maxBins - the number of bins to process from the raw data
        nTimeAvg - the number of raw rays to average in time when processing
        nRangeAvg - the number of bins to average in range when processing
        instrument - the name of the HSRL instrument providing raw data
        """
        QObject.__init__(self)

        # We do the hard work of product generation in separate 
        # ProductGeneratorThread threads.
        # Use exactly two threads since that provides a *slight* performance 
        # benefit over just one. More than two threads does not help because 
        # product generation is a compute-bound problem, and threading benefit 
        # under Python is limited in the compute-bound case by Python's 
        # Global Interpreter Lock (GIL), regardless of the number of processors
        # available.
        self.__threads = []
        for i in range(2):
            thread = ProductGeneratorThread(instrument, maxBins, nTimeAvg, 
                                            nRangeAvg)
            self.__threads.append(thread)
            # Re-emit newProductRay(ProductRay) signals from our worker threads
            # as newProductRay(ProductRay) signals from here
            QObject.connect(thread, SIGNAL("newProductRay(PyQt_PyObject)"), 
                            self, SIGNAL("newProductRay(PyQt_PyObject)"))
            thread.start()

        # We accumulate nTimeAvg raw rays before passing them to one of the
        # ProductGeneratorThread-s for processing.
        self.__nTimeAvg = nTimeAvg

        # list to accumulate raw rays before processing
        self.__rawRays = []

        self.__iwg1Socket = QUdpSocket()
        self.__iwg1Socket.bind(QHostAddress(QHostAddress.Any), 7071)
        QObject.connect(self.__iwg1Socket, SIGNAL("readyRead()"),
                        self.handleIWG1, Qt.QueuedConnection)

        self.__mtpSocket = QUdpSocket()
        self.__mtpSocket.bind(QHostAddress(QHostAddress.Any), 40002)
        QObject.connect(self.__mtpSocket, SIGNAL("readyRead()"),
                        self.handleMTP, Qt.QueuedConnection)

    def handleRawRay(self, rawRay):
        """
        Accept an incoming RawDataRay of raw accumulator data and process 
        computed products. Each time nTimeAvg raw rays are accumulated, 
        processing for a new ray of product data is initiated.
        """        
        # Append this ray to our raw ray list
        self.__rawRays.append(rawRay)

        # If we have accumulated enough raw data, pass our list of raw rays to
        # a ProductGeneratorThread to be turned into a ProductRay.
        if (len(self.__rawRays) == self.__nTimeAvg):
            for i in range(len(self.__threads)):
                # Pass processing on to the first thread which is not busy
                thread = self.__threads[i]
                if (not thread.isBusy()):
                    thread.generateProductRay(self.__rawRays)
                    break
                # If all threads are busy, drop the raw rays and complain
                if (i == (len(self.__threads) - 1)):
                    print "Product ray at", rawRay.rayTime(), \
                          "dropped because all threads are busy!"

            self.__rawRays = []


        return

    def handleIWG1(self):
        while self.__iwg1Socket.hasPendingDatagrams():
            size = self.__iwg1Socket.pendingDatagramSize()
            payload, hostaddr, val = self.__iwg1Socket.readDatagram(size)
            print 'got IWG1 packet'
        return;

    def handleMTP(self):
        while self.__mtpSocket.hasPendingDatagrams():
            size = self.__mtpSocket.pendingDatagramSize()
            payload, hostaddr, val = self.__mtpSocket.readDatagram(size)
            print 'got MTP packet'
        return;

from RawDataClient import RawDataClient
import signal

if __name__ == "__main__":
    def handleRaw(rawRay):
        print "main got a RawDataRay at", rawRay.rayTime()

    def handleProduct(pRay):
        print "main got a ProductRay at", pRay.rayTime(), "with products:", \
              pRay.productNames(), "in thread", thread.get_ident()
    # create a Qt core application
    print "Main thread is", thread.get_ident()
    app = QCoreApplication(sys.argv)

    # Client to read raw data
    reader = RawDataClient()

    # Product generator
    prodGen = HsrlProductGenerator(5, 4, 'gvhsrl')

    # Send raw data to the product generator
    QObject.connect(reader, SIGNAL("newRawDataRay(PyQt_PyObject)"), 
                    prodGen.handleRawRay, Qt.QueuedConnection)
    QObject.connect(reader, SIGNAL("newRawDataRay(PyQt_PyObject)"), 
                    handleRaw)

    # Handle products from the product generator with handleRay()
    QObject.connect(prodGen, SIGNAL("newProductRay(PyQt_PyObject)"), 
                    handleProduct, Qt.QueuedConnection)

    # Allow ^C to shut down the application
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    # start the raw data reader
    reader.start()

    # start the event loop
    sys.exit(app.exec_())
