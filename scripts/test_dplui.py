#!/usr/bin/env python

from PyQt4 import QtGui

from dplkit.ui.stream import Stream
from dplkit.ui.controller import GUIController
from dplkit.ui.widget import MplWidget
from dplkit.simple.filter import RollingWindowFilter,StructToDict,TransposeChannelsFilter

from hsrl.dpl_experimental.dpl_hsrl import dpl_hsrl
import hsrl.dpl_tools.substruct as substruct
import hsrl.dpl_tools.time_frame
from datetime import datetime,timedelta
import numpy
from matplotlib import colors
from matplotlib.ticker import LogFormatter

import sys

IMAGE_HISTORY_SIZE=240

def get_hsrl_frame_stream():
    starttime=datetime(2013,6,11,12,0,0)
    endtime=datetime(2013,6,11,14,0,0)
    #endtime=datetime(2013,6,12,0,0,0)

    # Get Librarian
    hsrllib=dpl_hsrl(instrument='bagohsrl')
    # Basic Narrator
    hsrlnar=hsrllib(start_time_datetime=starttime,end_time_datetime=endtime,min_alt_m=0.0,max_alt_m=20000.0, \
                timeres_timedelta=timedelta(seconds=30),altres_m=15) #calls search, returns a narrator
    altitudes = hsrlnar.altitudeAxis
    # Uncluster the complex frames
    hsrlnarsplitter=substruct.SubstructBrancher(hsrlnar) #returns another narrator that can extract nested streams

    hsrlinvnar=hsrlnarsplitter.narrateSubstruct('rs_inv') #inverted HSRL data nested stream

    # Get a per time-step frame stream
    hsrlinvframenar=hsrl.dpl_tools.time_frame.TimeGinsu(hsrlinvnar,'times')

    return StructToDict(hsrlinvframenar),altitudes

def get_hsrl_windowed_frame_stream():
    #starttime=datetime(2013,6,11,12,0,0)
    #endtime=datetime(2013,6,12,0,0,0)

    # Get Librarian
    hsrllib=dpl_hsrl(instrument='nshsrl', filetype="data")
    #hsrllib=dpl_hsrl(instrument='bagohsrl', filetype="data")
    # Basic Narrator
    hsrlnar=hsrllib(min_alt_m=0.0,max_alt_m=20000.0, reverse_padding=timedelta(seconds=45), \
                timeres_timedelta=timedelta(seconds=30),altres_m=15,window_width_timedelta=timedelta(seconds=60*60*2)) #calls search, returns a narrator
    #hsrlnar=hsrllib(start_time_datetime=starttime,end_time_datetime=endtime,min_alt_m=0.0,max_alt_m=20000.0, \
    #            timeres_timedelta=timedelta(seconds=30),altres_m=15) #calls search, returns a narrator
    altitudes = hsrlnar.altitudeAxis

    # Uncluster the complex frames
    hsrlnarsplitter=substruct.SubstructBrancher(hsrlnar) #returns another narrator that can extract nested streams

    hsrlinvnar=hsrlnarsplitter.narrateSubstruct('rs_inv') #inverted HSRL data nested stream

    # Get a per time-step frame stream
    hsrlinvframenar=hsrl.dpl_tools.time_frame.TimeGinsu(hsrlinvnar,'times')

    dict_stream = StructToDict(hsrlinvframenar)
    windowed_stream = RollingWindowFilter(dict_stream, window_size=IMAGE_HISTORY_SIZE)
    transposed_stream = TransposeChannelsFilter(windowed_stream)
    return transposed_stream,altitudes

def main():
    """Create a DPL UI and plot HSRL data
    """
    # Initialize the GUI application
    app = QtGui.QApplication([" "])

    # Create the HSRL stream
    print "Getting HSRL stream..."
    hsrl_frame_stream,altitudes = get_hsrl_frame_stream()
    print "Done"

    # Create the UI adapted stream (with artificial 1 second delay)
    ui_stream = Stream(hsrl_frame_stream, delay=1.0)

    # Create our plot widgets
    bs_line_plot = MplWidget(blit=True)
    bs_line_plot.setWindowTitle("HSRL Backscatter Perp.")
    fig = bs_line_plot.figure
    ax = fig.add_subplot(111)
    line = ax.plot(altitudes, numpy.array([numpy.nan]*1334))[0]
    ax.set_xlim(0, 20000)
    ax.set_ylim(-0.00002, 0.00002)
    ax.set_title("HSRL Perpendicular Backscatter (per-record)")
    ax.set_xlabel("Altitude (m)")
    ax.set_ylabel("Backscatter (1/(m sr))")
    ax.grid()

    # Create a controller to attach the frame stream to the widget
    controller = GUIController(ui_stream, bs_line_plot)
    controller.bind_line_y_to_channel("beta_a_backscat_perp", line)

    # Parallel Backscatter Plot
    bs_par_line_plot = MplWidget(blit=True)
    bs_par_line_plot.setWindowTitle("HSRL Backscatter Par.")
    fig = bs_par_line_plot.figure
    ax = fig.add_subplot(111)
    line = ax.plot(altitudes.squeeze(), numpy.array([numpy.nan]*1334))[0]
    ax.set_xlim(0, 20000)
    ax.set_ylim(-0.00002, 0.00002)
    ax.set_title("HSRL Parallel Backscatter (per-record)")
    ax.set_xlabel("Altitude (m)")
    ax.set_ylabel("Backscatter (1/(m sr))")
    ax.grid()

    # Create a controller to attach the frame stream to the widget
    bs_par_controller = GUIController(ui_stream, bs_par_line_plot)
    bs_par_controller.bind_line_y_to_channel("beta_a_backscat_par", line)

    # Show our window that we've made
    bs_line_plot.show()
    bs_par_line_plot.show()
    # Start the stream
    ui_stream.start()
    # Start the GUI application
    app.exec_()
    # Once the GUI is closed wait for our stream to end gracefully
    ui_stream.wait()

def main2():
    """Create a DPL UI and plot HSRL data
    """
    numpy.set_printoptions(linewidth=120)
    # Initialize the GUI application
    app = QtGui.QApplication([" "])

    # Create the HSRL stream
    print "Getting HSRL stream..."
    hsrl_frame_stream,altitudes = get_hsrl_windowed_frame_stream()
    print "Done"

    # Create the UI adapted stream (with artificial 1 second delay)
    ui_stream = Stream(hsrl_frame_stream, delay=0.5)
    #ui_stream = Stream(hsrl_frame_stream, delay=1.0)

    # Create our plot widgets
    bs_perp_image_plot = MplWidget(blit=True)
    bs_perp_image_plot.setWindowTitle("HSRL Backscatter Perp.")
    fig = bs_perp_image_plot.figure
    ax = fig.add_subplot(111)
    image = ax.imshow(numpy.zeros((1334,IMAGE_HISTORY_SIZE)), interpolation='bilinear', norm=colors.LogNorm(1e-8, 1e-3), aspect=IMAGE_HISTORY_SIZE/1334.0, animated=True)
    ax.set_title("HSRL Perpendicular Backscatter (per-record)")
    ax.set_xlabel("Records")
    ax.set_ylabel("Altitude (m)")
    ax.grid()
    log_f = LogFormatter(10, labelOnlyBase=False)
    fig.colorbar(image, ax=ax, format=log_f)

    # Create a controller to attach the frame stream to the widget
    controller = GUIController(ui_stream, bs_perp_image_plot)
    controller.bind_image_to_channel("beta_a_backscat_perp", image)

    # Parallel Backscatter Plot

    # Show our window that we've made
    bs_perp_image_plot.show()
    # Start the stream
    ui_stream.start()
    # Start the GUI application
    app.exec_()
    # Once the GUI is closed wait for our stream to end gracefully
    ui_stream.wait()

if __name__ == "__main__":
    sys.exit(main2())
    #sys.exit(main())
