#!/usr/bin/python
# -*- coding: utf-8 -*-
if __name__ == '__main__':

  from dpl_read_templatenetcdf import dpl_read_templatenetcdf
  import hsrl.data_stream.display_utilities as du
  import matplotlib.pyplot as plt
  import lg_base.core.open_config as oc
  import json

  plt.ion()
  v=dpl_read_templatenetcdf('outtest2.nc')
  
  instrument=v.raw_netcdf().instrument[:]

  (disp,conf)=du.get_display_defaults('dpl_plots.json','new')
  fd = oc.open_config('process_control.json')
  dd = json.load(fd)
  #self.rs_static.corr_adjusts = dd['corr_adjusts']
  process_defaults=dd['process_defaults']
  fd.close()

  for n in v:
    figs=du.show_images(instrument,n,None,{},[],process_defaults,disp,None,None,None)
    #print n.rs_inv.beta_a_backscat_par
    #print n.rs_inv.times
    # force a redraw
    firstfig=None
    plt.show(block=False)

    for x in figs:
        
        #      print 'updating  %d' % x.num
        
        
        fig = figs.figure(x)
        if firstfig==None or firstfig.number>fig.number:
            firstfig=fig
        
      # QApplication.processEvents()

        fig.set_size_inches(fig.get_size_inches(),forward=True)
        if hasattr(fig.canvas,"repaint"):
            fig.canvas.repaint()
        else:
            fig.canvas.draw()

    if firstfig!=None:
        firstfig.set_size_inches(firstfig.get_size_inches(),forward=True)
        if hasattr(firstfig.canvas,'repaint'):
            firstfig.canvas.repaint()
        else:
            firstfig.canvas.draw()

    raw_input("press a key")
