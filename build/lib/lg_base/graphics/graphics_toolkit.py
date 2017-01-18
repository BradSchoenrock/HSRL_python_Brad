#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.dates as mdate
from datetime import datetime,timedelta
import matplotlib.ticker as mticker
from pylab import get_cmap
import traceback
import sys,os
from matplotlib.font_manager import FontProperties
from matplotlib import colors 
import copy
from collections import OrderedDict,namedtuple
ClassColormapEntry=namedtuple('ClassColormapEntry','value name color')

def expLabels(base=10):
    def labelMaker(value,pos):
        return "$%i^{%i}$" % (base,value)
    def labelJSONMaker(value,pos):
        return "%i^%i" % (base,value)
    def labelExpMaker(value,pos):
        return "1e%i" % (value,)
    if 'tex' in os.getenv("GRAPHICSHACKS","").split(','):
        if base==10:
            return mticker.FuncFormatter(labelExpMaker)
        return mticker.FuncFormatter(labelJSONMaker)
    return mticker.FuncFormatter(labelMaker)

def logLabels(format,base=10):
    def labelMaker(value,pos):
        pval=base**value
        return format % pval
    return mticker.FuncFormatter(labelMaker)

def listLabels(liststrings,values=None):
    def labelMaker(value,pos):
        if values is None:
            v=int(value)
            if v<0:
                return "TOOLOW"
            if v>=len(liststrings):
                return "TOOHIGH"
            return liststrings[v]
        diff=[abs(value-x) for x in values]
        bestv=diff[0]
        besti=0
        for i,v in enumerate(diff):
            if v<bestv:
                bestv=v
                besti=i
        return liststrings[besti]
    if values is not None:
        assert(len(liststrings)==len(values))
    return mticker.FuncFormatter(labelMaker)

def verifyLogScale(display_defaults,plot_name,axis,minname=None,maxname=None):
    minv=display_defaults.get_value(plot_name,minname or (axis+' min'))
    if minv is not None and minv<=0.0:
        print 'ERROR: Logscale on '+axis+' axis explicitly set to negative min value on '+plot_name
        assert(minv is None or minv>0.0)
    maxv=display_defaults.get_value(plot_name,maxname or (axis+' max'))
    if maxv is not None and maxv<=0.0:
        print 'ERROR: Logscale on '+axis+' axis explicitly set to negative max value on '+plot_name
        assert(maxv is None or maxv>0.0)

class pyplot_figurelist:
    """ A container for creating and managing figures, with toolkit functions like redraw that may be used in a few places
    this version uses pyplot to create figures, which is useful if you use figure numbers
    """
    def __init__(self):
        self.figs=OrderedDict()
        import matplotlib.pyplot as plt
        self.plt=plt
        self.shown=set()

    def __iter__(self):
        return iter(self.figlist)

    def __contains__(self,v):
        return v in self.figs

    def __repr__(self):
        return 'figurelist'+repr(self.figlist)

    @property
    def figlist(self):
        return copy.copy(self.figs.keys())

    def showfigs(self,figlist):
        firstfig=None
        if not self.plt.isinteractive():
            self.plt.show(block=False)
        else:
            self.plt.show()

        for x in figlist:
            try:
                fig = self.figure(x)
            except:
                fig = x
                x=None
            #help(fig.canvas)
            try:
                self.redrawfig(fig,x)
            except ValueError:
                print 'Error with figure ',fig.number,x
            if firstfig is None or firstfig.number>fig.number:
                firstfig=fig
        if firstfig!=None:
            self.redrawfig(firstfig)

    def reorder(self,orderlist,remainingLast=False):
        oldfigs=self.figs
        self.figs=OrderedDict()
        if not remainingLast:
            for k,v in oldfigs.items():
                if k not in orderlist:
                    self.figs[k]=oldfigs.pop(k)
        for k in orderlist:
            if k in oldfigs:
                self.figs[k]=oldfigs.pop(k)
        if remainingLast:
            for k,v in oldfigs.items():
                self.figs[k]=v

    def shownew(self):
        l=[]
        for f in self:
            if f not in self.shown:
                l.append(f)
        self.showfigs(l)

    def showall(self):
        self.showfigs(self)

    def redrawfig(self,fig,figkey=None):
            #firstfig.set_size_inches(firstfig.get_size_inches(),forward=True)
            if not hasattr(fig,'canvas'):
                fig=self.figure(fig)
            if figkey is None:
                for k,f in self.figs.items():
                    if fig is f:
                        figkey=f
                        break                
            if figkey is not None:
                self.shown.add(figkey)
            i=fig.get_size_inches()*fig.get_dpi()
            fig.canvas.resize(i[0],i[1])
            fig.canvas.draw()
            if hasattr(fig.canvas,'repaint'):
                fig.canvas.repaint()


    def figure(self,val,*args,**kwargs):
        if not val in self:
            self.figs[val]=self.plt.figure(*args,**kwargs)
            self.figs[val].isNewFigure=True
        else:
           #plt.figure(self.figs[val].number)#still need to set to gcf... for now.
           #preference is to not imply a focus when multithreaded
           self.figs[val].isNewFigure=False
           f=self.figs[val]
           del self.figs[val]#reorder
           self.figs[val]=f
        self.shown.discard(val)
        return self.figs[val];

    def clear(self):
        for k in self:
            #help(self.figs[k])
            f=self.figs[k]
            f.clear()
            if hasattr(f,'fig_colorbars'):
                delattr(f,'fig_colorbars')
        self.shown.clear()
 
    def close(self,fig=None):
        if fig is not None:
            if not isinstance(fig,basestring):
                for k,f in self.figs.items():
                    if f is fig:
                        fig=k
                        break
            if fig not in self.figs:
                raise RuntimeError('Unknown figure '+repr(fig))
            self.figs[fig].clear()
            self.plt.close(self.figs[fig])
            self.shown.discard(fig)
            del self.figs[fig]
            return
        for k in self:
            self.close(k)
            #help(self.figs[k])
        #self.figs.clear()

    def __del__(self):
        self.close()

class matplotlib_figurelist: #for visual contexts, this is incomplete and technically useless. Agg only (for web)
    def __init__(self,forcedBackend=None):
        self.figs=OrderedDict()
        self.idx=0
        import matplotlib.figure as figplot
        self.figplot=figplot
        import matplotlib
        base=forcedBackend or matplotlib.get_backend()
        canvas=None
        if base=='Qt4Agg':
            import matplotlib.backends.backend_qt4agg as backend
            canvas=backend.FigureCanvasQTAgg
        elif base=='Qt4':
            import matplotlib.backends.backend_qt4 as backend
            canvas=backend.FigureCanvasQT
        elif base=='agg':
            import matplotlib.backends.backend_agg as backend
            canvas=backend.FigureCanvasAgg
        elif base=='MacOSX':
            import matplotlib.backends.backend_macosx as backend
            canvas=backend.FigureCanvasMac
        if canvas is None:
            raise NotImplementedError('Matplotlib Figurelist container using backend '+base)
        self.backendname=base
        self.backend=backend
        #self.show=self.backend.Show()
        self.canvas=canvas
        self.shown=set()

    def __iter__(self):
        return iter(self.figlist)

    def __contains__(self,v):
        return v in self.figs

    def __repr__(self):
        return 'figurelist'+repr(self.figlist)

    @property
    def figlist(self):
        return copy.copy(self.figs.keys())

    def reorder(self,orderlist,remainingLast=True):
        oldfigs=self.figs
        self.figs=OrderedDict()
        if not remainingLast:
            for k,v in oldfigs.items():
                if k not in orderlist:
                    self.figs[k]=oldfigs.pop(k)
        for k in orderlist:
            if k in oldfigs:
                self.figs[k]=oldfigs.pop(k)
        if remainingLast:
            for k,v in oldfigs.items():
                self.figs[k]=v

    def showfigs(self,figlist):
        firstfig=None

        for x in figlist:
            try:
                fig = self.figure(x)
            except:
                fig = x
                x=None
            #help(fig.canvas)
            self.redrawfig(fig,x)
            if firstfig is None or firstfig.number>fig.number:
                firstfig=fig
        if firstfig is not None:
            self.redrawfig(firstfig)

    def shownew(self):
        l=[]
        for f in self:
            if f not in self.shown:
                l.append(f)
        self.showfigs(l)

    def showall(self):
        self.showfigs(self)


    def redrawfig(self,fig,figkey=None):
            #firstfig.set_size_inches(firstfig.get_size_inches(),forward=True)
            if not hasattr(fig,'canvas'):
                fig=self.figure(fig)
            if figkey is None:
                for k,f in self.figs.items():
                    if fig is f:
                        figkey=f
                        break                
            if figkey is not None:
                self.shown.add(figkey)
            i=fig.get_size_inches()*fig.get_dpi()
            fig.canvas.resize(i[0],i[1])
            try:
                fig.canvas.draw()
                if hasattr(fig.canvas,'repaint'):
                    fig.canvas.repaint()
                try:
                    fig.canvas.flush_events()
                except NotImplementedError:
                    pass #nothing we can do about it
            except ValueError:
                print 'Error with figure ',fig.number,figkey

    def figure(self,val,*args,**kwargs):
        if not val in self:
            fig=self.figplot.Figure(*args,**kwargs)
            self.canvas(fig)
            fig.isNewFigure=True
            self.idx+=1
            fig.number=self.idx #fixme FIG NUMBERS are not safe
            fig.clear()
            self.figs[val]=fig
        else:
            #plt.figure(self.figs[val].number)#still need to set to gcf... for now.
            #preference is to not imply a focus when multithreaded
            self.figs[val].isNewFigure=False
            f=self.figs[val]
            del self.figs[val]
            self.figs[val]=f
        self.shown.discard(val)
        return self.figs[val];

    def clear(self):
        for k in self:
            f=self.figs[k]
            f.clear()
            if hasattr(f,'fig_colorbars'):
                delattr(f,'fig_colorbars')
        self.shown.clear()

    def close(self,fig=None):
        if fig is not None:
            if not isinstance(fig,basestring):
                for k,f in self.figs.items():
                    if f is fig:
                        fig=k
                        break
            if fig not in self.figs:
                raise RuntimeError('Unknown figure '+repr(fig))
            self.figs[fig].clear()
            del self.figs[fig]
            self.shown.discard(fig)
            return
        for k in self:
            self.close(k)

    def __del__(self):
        self.close()

figurelist=pyplot_figurelist

def plot_image(plot_name,instrument,times,image_name,color_map,title,defaults,figs,invert=False):
    """plot an image
       plot_image(plot_name
                 ,instrument
                 ,times
                 ,image_name
                 ,color_map
                 ,title
                 ,defaults
                 ,figs)"""

    f=figs.figure(plot_name,figsize=defaults.get_size('profile_graph_size'))
    f.canvas.set_window_title('Fig ' + str(f.number)
                + '        '+title)
    #f.set_size_inches(defaults.get_size('profile_graph_size'),forward=True)
    img = image_name.view(dtype=np.uint8)
    if invert:
        img=255-img
    if color_map is None:
        image = f.gca().imshow(img)
    else:
        image = f.gca().imshow(img,cmap = get_cmap(color_map))
    ax=f.gca()
    try:
        time_str = times[0].strftime('%d-%b-%Y %H:%M')
    except:
            try:
             time_str = times.strftime('%d-%b-%Y %H:%M')
            except:
             time_str = 'UNKNOWN TIME'
    title=instrument + ' ' + title + ' ' + time_str
    ax.set_title(title)
    image.set_clim(0, 255)
    
    return figs



def plot_mode_bits(plot_name,instrument,times,bits,bit_names,title
               ,display_defaults,figs):
    """plot_mode_bits(plot_name,instrument,times,bits,bit_names
           ,title,display_defaults,figs)
       plot mode bits as function of time
       plot_name         = plot name in display_defaults
       instrument        = instrument name, string appears in plot title
       times[ntimes]     = time vector used on x-axis
       bits[ntimes,nbits]= array of bits to display as function of time
       bit_names[nbits]  = names of bits
       title             = String that appears in plot title
       display_defaults  = structure from plot defaults *.json
       figs              = assigned by plot routines to number plots"""      


    verticalresizefactor=1.2 #used to prevent size change creep
    sz=list(display_defaults.get_size('image_size'))
    sz[1]*=verticalresizefactor or 1.0
    if 'figure' in display_defaults.get_labels(plot_name):
        fig_name=display_defaults.get_value(plot_name,'figure')#to override where this ends up being plotted
    else:
        fig_name=plot_name
    f=figs.figure(fig_name,dpi=int(display_defaults.get_value('image_pixels_per_inch','ppi',require=True)),figsize=sz)
    graph_setup(f,display_defaults,title,plot_name,True)
    bits=bits.copy()
    bits[bits==0] = np.NaN
    ax=f.gca()
    for i in range(bits.shape[1]):
        tmp=bits[:,i]
        tmp[np.isfinite(tmp)]=i
        if np.isfinite(tmp).any():
            ax.plot(times,tmp,linewidth=4)
            print 'plotting bit',i
        else:
            print 'skipping bit',i
    #for i in range(len(bit_names)):
    #f.gca().plot(times,reallist,linewidth=4)    
    ax.set_ylim((-0.5,len(bit_names)-0.5))
    make_xaxis(times, ax)

    ax.set_yticks(np.arange(0,len(bit_names)))
    ax.yaxis.set_major_formatter(listLabels(bit_names))# set_yticklabels(bit_names)
    ax.grid(True)  # ,linestyle='-')
    time_str = times[0].strftime('%d-%b-%Y')
    ax.set_title(instrument + ' '+title +'  '+time_str)
    #size=f.get_size_inches()
    #if verticalresizefactor!=None:
    #        f.set_size_inches([size[0],verticalresizefactor*size[1]])

def plot_vs_altitude(plot_name,instrument,times,altitudes,x_vars,colors,widths
                     ,legend_list,legend_position,xlabel,xunits,title,clear_fig
                     ,display_defaults,figs,yunits=None):
        """plot variables vs alt with altitudes defined by a common alt vector
         plot_vs_altitude(name     #plot name contained in plot_defaults *.json
                 ,instrument       #instrument name,string appears in plot title
                 ,time             #matplotlib time this appears in plot title
                 ,altitudes        #altitude vector (meters)
                 ,[var 1,var2, ..] #vars to plot, specified at alt vector pts
                 ,['r','b',....]   #colors, None default colors
                 ,[3, 1]           #widths, None sets widths = 2
                 ,['p1','p2']      #legend list, None = no legend
                 ,'upper left'     #legend position, None ok if legend list []
                 ,'Pressure'       #xlabel
                 ,'mb'             #x units
                 ,'etalon pressure'#string that appears in plot title
                 ,clear_fig        #=1, clr fig before ploting, else plot over
                 ,display_defaults #structure from plot defaults *.json         
                 ,figs)"""
               
                       
        verticalresizefactor=1.2 #used to prevent size change creep
        if yunits is None:
            yunits='km'
            altitudes=altitudes/1000.0
        sz=list(display_defaults.get_size('profile_graph_size'))
        sz[1]*=verticalresizefactor or 1.0
        f=figs.figure(plot_name,figsize=sz)
        if clear_fig:
                f.clf()
        if f.isNewFigure:
            f.canvas.set_window_title('Fig ' + str(f.number)
                + '       ' +title)
            #f.set_size_inches(display_defaults.get_size('profile_graph_size'),forward=True)
       
            f.subplots_adjust(
                top=.9,
                bottom=.1,
                left=.1,
                right=.9,
                hspace=0,
                wspace=0,
                )
            
        if colors is not None and len(colors)==0:
            print 'COlors should be None'
            colors = None
        if colors is None:
            colors=['r','b','g','c','k','m','y','k','r','b','g','c','k','m','y']
          
        if widths is not None and len(widths)==0:
            print 'widths should be None'
            widths = None
        if xunits is not None and len(xunits)==0:
            print 'xunits should be None'
            xunits = None
        if legend_list is not None and len(legend_list)==0:
            print 'legend_list should be None'
            legend_list = None
        if widths is None:
            widths= [2]*len(x_vars)
        try:
            for i in range(len(x_vars)):
                if not display_defaults.get_value(plot_name,'log_linear') is None \
                       and display_defaults.get_value(plot_name,'log_linear').find('log')>=0:
                    #eliminate values <=0 for log plot
                    temp = x_vars[i].copy()
                    temp[np.isnan(temp)]=0.0
                    temp[temp<=0.0 ] = np.NaN
                    f.gca().plot(temp,altitudes,color=colors[i],linewidth=widths[i])
                else:
                    f.gca().plot(x_vars[i],altitudes,color=colors[i],linewidth=widths[i])
        except:
            print ' '
            print '****plot_vs_altitude was unable to create "'+plot_name+ '"   **************'
            print 'len(altitudes) = ',len(altitudes)
            print 'x_var[0].shape = ',x_vars[0].shape
            print 'len(x_var)= ',len(x_vars)  
            print 'colors' ,colors
            print 'widths', widths
            traceback.print_exc()
            return
        
        ax=f.gca()
        [xlow,xhigh]=ax.get_xlim()
        if display_defaults.get_value(plot_name,'x min') is not None:    
            ax.set_xlim(xmin=(display_defaults.get_value(plot_name,'x min')))
        if display_defaults.get_value(plot_name,'x max') is not None:
            ax.set_xlim(xmax=(display_defaults.get_value(plot_name,'x max')))
        [xlow,xhigh]=ax.get_xlim()
        if xunits is not None:
            ax.set_xlabel(xlabel+' ('+xunits+')')
        else:
            ax.set_xlabel(xlabel)
       
        if not display_defaults.get_value(plot_name,'log_linear') is None:
            #if display_defaults.get_value(plot_name,'log_linear').find('log')>=0:
            if 'log' in display_defaults.get_value(plot_name,'log_linear'):
                try:
                    verifyLogScale(display_defaults,plot_name,'x')
                    ax.set_xscale('log')
                except ValueError:
                    ax.set_xscale('linear')
                    print "Can't set xscale to log"
                    pass
        ax.set_ylabel('MSL Altitudes ('+yunits+')')
        ax.yaxis.set_units(yunits)
        ax.grid(True)
        try:
            if times[0].strftime('%d') == times[-1].strftime('%d'):
                time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                       +' --> '+ times[-1].strftime('%H:%M')
            else:
                time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                      + ' --> ' + times[-1].strftime('%d-%b %H:%M')
        except:
            try:
             time_str = times.strftime('%d-%b-%Y %H:%M')
            except:
             time_str = 'UNKNOWN TIME'
            
        ax.set_title(instrument + ' '+title +'  '
                  + time_str)
        if legend_list is not None:
           fontP = FontProperties()
           fontP.set_size('small')
           ax.legend(legend_list,loc=legend_position,prop = fontP)
           #ax.legend(legend_list,legend_position)
        #size=f.get_size_inches()
        #if verticalresizefactor!=None:
        #    f.set_size_inches([size[0],verticalresizefactor*size[1]])
        return figs

def plot_map_latlon(plot_name,instrument,times,lon_vars,lat_vars,colors,markers,marker_size
            ,line_style,widths,title,text_str,text_position_lon,text_position_lat,text_angle,display_defaults,figs,fig_name=None):
    """create an x-y plot with discreate points and/or lines including legends and arbitrary text
    
       plot_xy(plot_name      #plot name as contained in plot_defaults *.json
           ,instrument        #instrument name,  this string appears in plot title
           ,time              #matplot vector time, this time appears in plot title
           ,x_vars            #list of x_vectors to plot    e.g. [[3,5],[6,7]...]
           ,y_vars            #list of y_vectors to plot    e.g. [[1,2],[3,4]...]
           ,colors            #list of colors for points,   e.g. ['k','r'.......]
           ,markers           #list of markers for points,  e.g. ['*', 'o' .....]
           ,marker_size       #list of marker sizes         e.g. [1, 2 .........]
           ,line_style        #list of linestyles,          e.g. ['_','o'.......]
           ,widths            #list of linewidths           e.g. [None,1,2...........]
           ,legend_list       #legend list (one for each variable)
           ,legend_position   # e.g. 'upper left'
           ,xlabel            #x-axis label (string) 
           ,x_units           #x units string
           ,ylabel            #y-axis label (string)
           ,y_units           #y units string
           ,title             #a string that will appear in plot title
           ,text_str          #list of strings to place on plot
           ,text_position_x   #list of x-positions for text_str entries
           ,text_position_y   #list of y-positions for text_str entries
           ,text_angle        #list of text angles for text_str entries
           ,display_defaults  #dictonary read from display_defaults.json
           ,figs)""" 
    
    if markers is not None and len(markers)==0:
        print 'MARKERS SHOULD BE NONE'
        markers=None
    if markers is None:
        markers=['None']*len(lat_vars)
        marker_size=[0]*len(lat_vars)
    elif len(markers) == 1:
        markers = markers*len(lat_vars)
        marker_size=[marker_size[0]]*len(lat_vars)
    if colors is not None and len(colors)==0:
        print 'colors SHOULD BE NONE'
        colors=None
    if text_str is not None and len(text_str)==0:
        print 'text_str SHOULD BE NONE'
        text_str=None
    if widths is not None and len(widths)==0:
        print 'widths SHOULD BE NONE'
        widths=None

    if colors is None:
        colors=['r','b','g','c','k','m','y','k','r','b','g','c','k','m','y']
        if len(colors)>15:
            colors=colors+['k']*(len(colors)-14)
    elif len(colors) == 1 :
        colors = colors*(len(lat_vars))
    if len(line_style) == 1:
        line_style = line_style*len(lat_vars)
    if len(widths) ==1:
        widths = widths*len(lat_vars)
    n_vars =len(lat_vars)
    t_vars = 0 if text_str is None else len(text_str)
    if widths is None:
        widths=[]
        for i in range(n_vars):
            widths.append(0.0)
        
    subplotval=0
    if fig_name is not None:
        pass
    elif 'figure' in display_defaults.get_labels(plot_name):
        fig_name=display_defaults.get_value(plot_name,'figure')#to override where this ends up being plotted
    else:
        fig_name=plot_name
    aspect=display_defaults.get_size('profile_graph_size')
   
    try:
        aspect=float(aspect[1])/aspect[0]
    except:
        aspect=1.0
   
    f=figs.figure(fig_name,figsize=display_defaults.get_size('profile_graph_size'))
    if 'subplot' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
        subplotval=display_defaults.get_value(plot_name,'subplot')
    
    if subplotval>0:
        if 'rect' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
            rect=display_defaults.get_value(plot_name,'rect')
            f.add_axes(rect)
        else:
            f.add_subplot(subplotval)
    else:
        ax=[.08,.05,0,0]#,.96,.95]
        ax[2]=1.0-(2*ax[0])
        ax[3]=1.0-(2*ax[1])-.05#vertical room for text
        f.add_axes(ax)
        f.canvas.set_window_title('Fig ' + str(f.number)
                + '        '+title)
    #f.set_size_inches(display_defaults.get_size('profile_graph_size'),forward=True)
    from mpl_toolkits.basemap import Basemap
    lonmin=None
    lonmax=None
    latmin=None
    latmax=None
   
    for i in range(len(lat_vars)):
        _lonmin=float(np.nanmin(lon_vars[i]))
        _lonmax=float(np.nanmax(lon_vars[i]))
        _latmin=float(np.nanmin(lat_vars[i]))
        _latmax=float(np.nanmax(lat_vars[i]))
        lonmin=_lonmin if lonmin is None else np.nanmin([_lonmin,lonmin])
        lonmax=_lonmax if lonmax is None else np.nanmax([_lonmax,lonmax])
        latmin=_latmin if latmin is None else np.nanmin([_latmin,latmin])
        latmax=_latmax if latmax is None else np.nanmax([_latmax,latmax])

    latpad=max([.1,(latmax-latmin)/2.0])
    latmax+=latpad
    latmin-=latpad
    lonpad=max([.1,(lonmax-lonmin)/2.0])
    lonmax+=lonpad
    lonmin-=lonpad
    casp=(latmax-latmin)/(lonmax-lonmin)
    if casp<aspect:
        latheight=aspect*(lonmax-lonmin)
        latmid=(latmax-latmin)/2.0+latmin
        latmax=latmid+latheight/2.0
        latmin=latmid-latheight/2.0
    else:
        lonheight=(1.0/aspect)*(latmax-latmin)
        lonmid=(lonmax-lonmin)/2.0+lonmin
        lonmax=lonmid+lonheight/2.0
        lonmin=lonmid-lonheight/2.0
  
    
    if display_defaults.get_value(plot_name,'x min') is not None:
        lonmin=display_defaults.get_value(plot_name,'x min')
    if display_defaults.get_value(plot_name,'x max') is not None:
        lonmax=display_defaults.get_value(plot_name,'x max')
    if display_defaults.get_value(plot_name,'y min') is not None:
        latmin=display_defaults.get_value(plot_name,'y min')
    if display_defaults.get_value(plot_name,'y max') is not None:
        latmax=display_defaults.get_value(plot_name,'y max')
    if latmax >90.0 or latmin <-90.0:
        print
        print
        print
        print 'latitude out of range--probably missing latitude array when expecting mobile platform'
        print 'setting latmax = 90.0, latmin = -90.0'
        print
        print
        latmax =90.0
        latmin = -90.0
    if lonmax >720 or lonmin <-360.0:
        print
        print
        print
        print 'longitude out of range--probably missing longitude array when expecting mobile platform'
        print 'setting lonmax = 360.0, lonmin = 0.0'
        print
        print
        lonmax =360.0
        lonmin = 0.0
        
    latheight=(latmax-latmin)
    if latheight<.125:
        res=.01
    elif latheight<1:
        res=.1
    elif latheight<2.5:
        res=.25
    elif latheight<5:
        res=1
    elif latheight<10:
        res=2.5
    else:
        res=5
       
    m = Basemap(lonmin,latmin,lonmax,latmax,ax=f.gca(),resolution='i',area_thresh=.1)#,projection='merc')
    try:
        for i in range(len(lat_vars)):
           m.plot(lon_vars[i],lat_vars[i],color=colors[i%len(colors)],marker=markers[i%len(markers)]
                         ,markersize=marker_size[i%len(marker_size)],linestyle=line_style[i%len(line_style)]
                         ,linewidth=float(widths[i%len(widths)]),latlon=True)
        #formatter = mticker.ScalarFormatter(useOffset=False)
        #f.gca().yaxis.set_major_formatter(formatter)
        #f.gca().xaxis.set_major_formatter(formatter)
                   
    except Exception, err:
            print ' '
            print '****plot_vs_xy was unable to create "'+plot_name+ '"   **************' 
            sys.stderr.write('ERROR: %s\n' % str(err))

            for i in range(len(lat_vars)):
                print 'i=',i
                print 'x_var[i].shape = ',lon_vars[i].shape
                print 'y_var[i].shape = ',lat_vars[i].shape
                print 'colors[%i]' % (i%len(colors)) ,colors[i%len(colors)]
                print 'linestyle[%i]' % (i%(len(line_style))),line_style[i%len(line_style)]
                print 'line_widths[%i]' % (i%(len(widths))),widths[i%len(widths)]
                print 'markers[%i]' % (i%(len(markers))), markers[i%len(markers)]
                print 'marker_size[%i]' % (i%(len(marker_size))),marker_size[i%len(marker_size)]
#                print 'text_str',text_str[i]
#                print 'text_position_x',text_position_x[i]
#                print 'text_position_y',text_position_y[i]
            traceback.print_exc()
            return

    #m.drawlsmask()
    m.drawcoastlines()
    m.drawrivers()
    m.fillcontinents()
    m.drawstates()
    m.drawcountries()
    #for sc in [1,2,5,10,25]:
    #    try:
    #        m.warpimage(scale=sc)
    #        break
    #    except:
    #        pass
    m.drawmeridians(np.arange(0.0,360.0,res),labels=[True,False,False,True],fmt='%g')
    m.drawparallels(np.arange(-90.0,90.0,res),labels=[False,True,True,False],fmt='%g')
   
    ax=f.gca()    

    try:
        if   times[0].strftime('%d') == times[-1].strftime('%d'): 
            time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                     +' --> ' + times[-1].strftime('%H:%M')
        else:
            time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                     +' --> ' + times[-1].strftime('%d-%b %H:%M')
    except: 
            try:
             time_str = times.strftime('%d-%b-%Y %H:%M')
            except:
             time_str = 'UNKNOWN TIME'
    ax.set_title(instrument + ' '+title +'  '
                  + time_str)
    
    if text_str is not None:
        text_x,text_y=m(text_position_lon,text_position_lat)
        for i in range(len(text_str)):
             if 0:
                 print ' '
                 print plot_name
                 print 'text=',text_str[i]
                 print 'lon=',text_position_lon[i]
                 print 'lat=',text_position_lat[i]
                 print 'x=',text_x[i]
                 print 'y=',text_y[i]
                 print 'theta=',text_angle[i]
             ax.text(
                    text_x[i],
                    text_y[i],
                    text_str[i],
                    fontsize=14,
                    rotation=text_angle[i],
                    verticalalignment='center',
                    horizontalalignment='center',
                    )
    return figs    

def unique_mask(_x,_y):
    hasset=set()
    x=_x.ravel()
    y=_y.ravel()
    #print 'unique',x.size
    ret=np.logical_and(np.isfinite(x),np.isfinite(y))
    #print x.shape
    for idx in range(x.size):
        if not ret[idx].all():
            continue
        st='%g %g' % (float(x[idx]),float(y[idx]))
        #print st
        if st not in hasset:
            hasset.add(st)
            continue
        ret[idx]=False
    return ret.reshape(_x.shape)

def plot_xy(plot_name,instrument,times,x_vars,y_vars,colors,markers,marker_size
            ,line_style,widths,legend_list,legend_position,xlabel,xunits
            ,ylabel,yunits,title,text_str,text_position_x,text_position_y,text_angle,display_defaults,figs,fig_name=None):
    """create an x-y plot with discreate points and/or lines including legends and arbitrary text
    
       plot_xy(plot_name      #plot name as contained in plot_defaults *.json
           ,instrument        #instrument name,  this string appears in plot title
           ,time              #matplot vector time, this time appears in plot title
           ,x_vars            #list of x_vectors to plot    e.g. [[3,5],[6,7].....]
           ,y_vars            #list of y_vectors to plot    e.g. [[1,2],[3,4].....]
           ,colors            #list of colors for points,   e.g. ['k','r'.........]
           ,markers           #list of markers for points,  e.g. ['*', 'o' .......]
           ,marker_size       #list of marker sizes         e.g. [1, 2 ...........]
           ,line_style        #list of linestyles,          e.g. ['None','_','o'..]
           ,widths            #list of linewidths           e.g. [1,2.............]
           ,legend_list       #legend list (one for each variable)
           ,legend_position   # e.g. 'upper left'
           ,xlabel            #x-axis label (string) 
           ,x_units           #x units string
           ,ylabel            #y-axis label (string)
           ,y_units           #y units string
           ,title             #a string that will appear in plot title
           ,text_str          #list of strings to place on plot
           ,text_position_x   #list of x-positions for text_str entries
           ,text_position_y   #list of y-positions for text_str entries
           ,text_angle        #list of text angles for text_str entries
           ,display_defaults  #dictonary read from display_defaults.json
           ,figs)"""
 
    if widths is not None and len(widths)==0:
        print 'widths SHOULD BE NONE'
        widths=None
    if markers is not None and len(markers)==0:
        print 'markers SHOULD BE NONE'
        markers=None
    if colors is not None and len(colors)==0:
        print 'colors SHOULD BE NONE'
        colors=None
    if text_str is not None and len(text_str)==0:
        print 'text_str SHOULD BE NONE'
        text_str=None
    if legend_list is not None and len(legend_list)==0:
        print 'legend_list SHOULD BE NONE'
        legend_list=None
    
    if markers is None:
        markers=['None']*len(y_vars)
        marker_size=[0]*len(y_vars)
    elif len(markers) == 1:
        markers = markers*len(y_vars)
        marker_size=[marker_size[0]]*len(y_vars)
    if colors is None:
        colors=['r','b','g','c','k','m','y','k','r','b','g','c','k','m','y']
        if len(colors)>15:
            colors=colors+['k']*(len(colors)-14)
    elif len(colors) == 1 :
        colors = colors*(len(y_vars))
    if len(line_style) == 1:
        line_style = line_style*len(y_vars)
    if len(widths) ==1:
        widths = widths*len(y_vars)
    n_vars =len(y_vars)
    t_vars = 0 if text_str is None else len(text_str)
    if widths is None:
        for i in range(n_vars):
            widths.append(0.0)
        
    subplotval=0
    if fig_name is not None:
        pass
    elif 'figure' in display_defaults.get_labels(plot_name):
        fig_name=display_defaults.get_value(plot_name,'figure')#to override where this ends up being plotted
    else:
        fig_name=plot_name
    f=figs.figure(fig_name,figsize=display_defaults.get_size('image_size'))
    if 'subplot' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
        subplotval=display_defaults.get_value(plot_name,'subplot')
    
    if subplotval>0:
        if 'rect' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
            rect=display_defaults.get_value(plot_name,'rect')
            f.add_axes(rect)
        else:
            f.add_subplot(subplotval)
    else:
        pass
        """
        f.add_axes([.08,.05,.96,.95])
        """
        f.canvas.set_window_title('Fig ' + str(f.number)
                + '        '+title)
               
    #f.set_size_inches(display_defaults.get_size('profile_graph_size'),forward=True)
   
    try:
        for i in range(len(y_vars)):
            if line_style[i%len(line_style)]=='None' or float(widths[i%len(widths)])==0:
                m=unique_mask(x_vars[i],y_vars[i])
                x_vars[i]=x_vars[i][m]
                y_vars[i]=y_vars[i][m]
            f.gca().plot(x_vars[i],y_vars[i],color=colors[i%len(colors)],marker=markers[i%len(markers)]
                         ,markersize=marker_size[i%len(marker_size)],linestyle=line_style[i%len(line_style)]
                         ,linewidth=float(widths[i%len(widths)]) )
        formatter = mticker.ScalarFormatter(useOffset=False)
        f.gca().yaxis.set_major_formatter(formatter)
        f.gca().xaxis.set_major_formatter(formatter)
               
    except Exception, err:
            traceback.print_exc()
            print ' '
            print '****plot_vs_xy was unable to create "'+plot_name+ '"   **************' 
            sys.stderr.write('ERROR: %s\n' % str(err))

            for i in range(len(y_vars)):
                print 'i=',i
                print 'x_var[i].shape = ',x_vars[i].shape
                print 'y_var[i].shape = ',y_vars[i].shape
                print 'colors[%i]' % (i%len(colors)) ,colors[i%len(colors)]
                print 'linestyle[%i]' % (i%(len(line_style))),line_style[i%len(line_style)]
                print 'line_widths[%i]' % (i%(len(widths))),widths[i%len(widths)]
                print 'markers[%i]' % (i%(len(markers))), markers[i%len(markers)]
                print 'marker_size[%i]' % (i%(len(marker_size))),marker_size[i%len(marker_size)]
#                print 'text_str',text_str[i]
#                print 'text_position_x',text_position_x[i]
#                print 'text_position_y',text_position_y[i]
            traceback.print_exc()
            return
    ax=f.gca()    
    #[xlow,xhigh]=ax.get_xlim()
    #[ylow,yhigh]=ax.get_ylim()
    if not display_defaults.get_value(plot_name,'y max') is None:
            yhigh = display_defaults.get_value(plot_name,'y max')
            ax.set_ylim(ymax=yhigh)
    if not display_defaults.get_value(plot_name,'y min') is None:
            ylow = display_defaults.get_value(plot_name,'y min')
            ax.set_ylim(ymin=ylow)
    if not display_defaults.get_value(plot_name,'x max') is None:
        xhigh = display_defaults.get_value(plot_name,'x max')
        ax.set_xlim(xmax=xhigh)
    if not display_defaults.get_value(plot_name,'x min') is None:
        xlow = display_defaults.get_value(plot_name,'x min')  
        ax.set_xlim(xmin=xlow)
    #ax.set_xlim([xlow,xhigh])
    #ax.set_ylim([ylow, yhigh])
    if display_defaults.get_value(plot_name,'x log'):
            try:
                verifyLogScale(display_defaults,plot_name,'x')
                ax.set_xscale('log')
            except ValueError:
                ax.set_xscale('linear')
                print "Warning: could not use log scale for ", title
    if display_defaults.get_value(plot_name,'y log'):
            try:
                verifyLogScale(display_defaults,plot_name,'y')
                ax.set_yscale('log')
            except ValueError:
                print "Warning: could not use log scale for ", title
    if yunits is not None:
        ax.set_ylabel(ylabel+' ('+yunits+')')
        ax.yaxis.set_units(yunits)
    else:
        ax.set_ylabel(ylabel)
    if xunits is not None:    
        ax.set_xlabel(xlabel+' ('+xunits+')')
        ax.xaxis.set_units(xunits)
    else:
        ax.set_xlabel(xlabel)
    ax.grid(True)
    try:
        if   times[0].strftime('%d') == times[-1].strftime('%d'): 
            time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                     +' --. ' + times[-1].strftime('%H:%M')
        else:
            time_str = times[0].strftime('%d-%b-%Y %H:%M')\
                     +' --> ' + times[-1].strftime('%d-%b %H:%M')
    except:
            try:
             time_str = times.strftime('%d-%b-%Y %H:%M')
            except:
             time_str = 'UNKNOWN TIME'
    ax.set_title(instrument + ' '+title +'  '
                  + time_str)
    if legend_list is not None:
        fontP = FontProperties()
        fontP.set_size('small')
        ax.legend(legend_list,loc=legend_position,prop = fontP)
        #ax.legend(legend_list,legend_position)
    
    if text_str is not None:
        for i in range(len(text_str)):
             if 0:
                 print ' '
                 print plot_name
                 print 'text=',text_str[i]
                 print 'x=',text_position_x[i]
                 print 'y=',text_position_y[i]
                 print 'theta=',text_angle[i]
             ax.text(
                    text_position_x[i],
                    text_position_y[i],
                    text_str[i],
                    fontsize=14,
                    rotation=text_angle[i],
                    verticalalignment='center',
                    horizontalalignment='center',
                    )     
    return figs    


def plot_vs_x(plot_name,instrument,times,x_vars,y_vars,colors,widths,legend_list \
                ,legend_position,xlabel,xunits,ylabel,yunits,title,clear_fig,display_defaults,figs):
    """plot one or more vairables vs x with variables  specified
        by a common x vector.
        
        plot_vs_x(name                    #plot name as contained in plot_defaults *.json
                 ,instrument                 #instrument name,string that appears in plot title
                 ,x_vector                   #x vector
                 ,[variable 1,variable2, ..] #list of variables to plot, specified at time vector pts
                 ,['r','b',....]             #list of colors, [] default colors
                 ,[3, 1]                     #list of widths, [] sets widths = 2
                 ,['p1','p2']                #legend list, [] = no legend
                 ,'upper left'               #legend position, [] ok if legend list []
                 ,'Current '                 #xlabel
                 ,'A'                        #x units
                 ,'pressure'                 #ylabel
                 ,'mb'                       #y units
                 ,'etalon pressure'          #string that appears in plot title
                 ,False                      # if true, will clear and setup figure
                 ,display_defaults           # structure obtained from plot defaults *.json         
                 ,figs)"""                
       
    if colors is not None and len(colors)==0:
        print 'colors SHOULD BE NONE'
        colors=None
    if widths is not None and len(widths)==0:
        print 'widths SHOULD BE NONE'
        widths=None
    if yunits is not None and len(yunits)==0:
        print 'yunits SHOULD BE NONE'
        yunits=None
    if xunits is not None and len(xunits)==0:
        print 'xunits SHOULD BE NONE'
        xunits=None
    if legend_list is not None and len(legend_list)==0:
        print 'legend_list SHOULD BE NONE'
        legend_list=None
    f=figs.figure(plot_name,dpi=int(display_defaults.get_value('image_pixels_per_inch','ppi')),figsize=display_defaults.get_size('image_size'))
    graph_setup(f, display_defaults,title,plot_name,clear_fig)
    if colors==[]:
        colors=['r','b','g','c','k','m','y','k','r','b','g','c','k','m','y','r','b'\
                ,'g','c','k','m','y']
          
    if widths == []:
        widths= [2]*len(y_vars)
    

    for i in range(len(y_vars)):
        if widths[i] == 0:
            marker = 'o'
        else:
            marker = None

    try:        
            f.gca().plot(x_vars,y_vars[i],color=colors[i],linewidth=widths[i],marker=marker)
    except:
        print ' '
        print '****plot_vs_x was unable to create "'+plot_name+ '"   **************'
        print 'x_vars[0].shape = ',x_vars[0].shape
        print 'y_var[0].shape = ',y_vars[0].shape
        print 'len(y_var)= ',len(y_vars)  
        print 'colors' ,colors
        print 'widths', widths
        traceback.print_exc()
        return
        
    ax=f.gca()
   
       

    
    #if provided in display_defaults, reset xlimits 
    [xlow,xhigh]=ax.get_xlim()
    [ylow,yhigh]=ax.get_ylim()
    if not display_defaults.get_value(plot_name,'x max') is None:
            xhigh = np.float(display_defaults.get_value(plot_name,'x max'))
            ax.set_xlim(right=xhigh)
    if not display_defaults.get_value(plot_name,'x min') is None:
            xlow = np.float(display_defaults.get_value(plot_name,'x min'))
            ax.set_xlim(left=xlow)
    if not display_defaults.get_value(plot_name,'high') is None:
        yhigh = display_defaults.get_value(plot_name,'high')
        ax.set_ylim([ylow, yhigh])
    if not display_defaults.get_value(plot_name,'low') is None:
        ylow = display_defaults.get_value(plot_name,'low')  
        ax.set_ylim(ylow,yhigh)

    if not display_defaults.get_value(plot_name,'log_linear') is None \
              and display_defaults.get_value(plot_name,'log_linear').find('log')>=0:
        try:
            verifyLogScale(display_defaults,plot_name,'y','low','high')
            ax.set_yscale('log')
        except ValueError:
            print "Warning: could not use log scale for ", title
    
   
    if yunits is not None:
        ax.set_ylabel(ylabel+' (' + yunits + ')')
        ax.yaxis.set_units(yuints)
    else:
        ax.set_ylabel(ylabel)
    if xunits is not None:
        ax.set_xlabel(xlabel+' ('+xunits + ')')
        ax.xaxis.set_units(xunits)
    else:    
        ax.set_xlabel(xlabel)
    time_str = times[0].strftime('%d-%b-%Y %H:%M') 
    ax.set_title(instrument + ' '+title +'  '
                  + time_str)
    if legend_list is not None:
        ax.legend(legend_list,loc=legend_position)
    ax.grid(True)
    return figs
        
def plot_vs_time(plot_name,instrument,times,y_vars,colors,widths,legend_list\
                ,legend_position,ylabel,yunits,title,clear_fig,display_defaults,figs):
    print 'rendering ',plot_name
    """plot one or more vairables vs time with variables specified at times specified
        by a common time python datetime vector.
        
        plot_vs_time(name                    #plot name as contained in plot_defaults *.json
                 ,instrument                 #instrument name,string that appears in plot title
                 ,times                      #python datetime vector
                 ,[variable 1,variable2, ..] #list of variables to plot, specified at time vector pts
                 ,['r','b',....]             #list of colors, [] default colors
                 ,[3, 1]                     #list of widths, [] sets widths = 2, if 0 no lines "o"
                 ,['p1','p2']                #legend list, [] = no legend
                 ,'upper left'               #legend position, [] ok if legend list []
                 ,'Pressure'                 #ylabel
                 ,'mb'                       #y units
                 ,'etalon pressure'          #string that appears in plot title
                 ,False                      # if true, will clear and rebuild fig
                 ,display_defaults           # structure obtained from plot defaults *.json         
                 ,figs)"""                
       
    if colors is not None and len(colors)==0:
        print 'colors SHOULD BE NONE'
        colors=None
    if widths is not None and len(widths)==0:
        print 'widths SHOULD BE NONE'
        widths=None
    if yunits is not None and len(yunits)==0:
        print 'yunits SHOULD BE NONE'
        yunits=None
    if legend_list is not None and len(legend_list)==0:
        print 'legend_list SHOULD BE NONE'
        legend_list=None

    needXLabel=True
    if 'figure' in display_defaults.get_labels(plot_name):
        fig_name=display_defaults.get_value(plot_name,'figure')#to override where this ends up being plotted
    else:
        fig_name=plot_name
    f=figs.figure(fig_name,dpi=int(display_defaults.get_value('image_pixels_per_inch','ppi',require=True)),figsize=display_defaults.get_size('image_size'))

    graph_setup(f, display_defaults,title,plot_name,clear_fig)
    if colors is None:
        colors=['r','b','g','c','k','m','y','k']
    if widths is None:
        widths= [2]*len(y_vars)
    log_ok = True
    try:
        for i in range(len(y_vars)):
            if widths[i%len(widths)] == 0:
                marker = 'o'
            else:
                marker = None
            if display_defaults.get_value(plot_name,'log_linear')=='log':
                 temp = y_vars[i].copy()
                 temp[temp <= 0] = np.NaN
                 f.gca().plot(times,temp,color=colors[i % len(colors)],linewidth=widths[i % len(widths)],marker=marker)
             
            else:
                f.gca().plot(times,y_vars[i],color=colors[i % len(colors)],linewidth=widths[i % len(widths)],marker=marker)
           
    except:
        print ' '
        print '****plot_vs_time was unable to create "'+plot_name+ '"   **************'
        print 'len(times) = ',len(times)
        print 'y_var[0].shape = ',y_vars[0].shape
        print 'len(y_var)= ',len(y_vars)  
        print 'colors' ,colors
        print 'widths', widths
        print 'y_vars[%i]'%i
        print y_vars[i]
        traceback.print_exc()
        return
        
    ax=f.gca()

    make_xaxis(times, ax)
    if not display_defaults.get_value(plot_name,'log_linear') is None \
          and display_defaults.get_value(plot_name,'log_linear') =='log': 
        try:
            verifyLogScale(display_defaults,plot_name,'y')
            ax.set_yscale('log')    
        except ValueError:
            print "Warning: could not use log scale for ", title
            ax.set_yscale('linear')    
    [ylow,yhigh]=ax.get_ylim()
    doSetLim=False
    if not display_defaults.get_value(plot_name,'y max') is None:
        yhigh = display_defaults.get_value(plot_name,'y max')
        doSetLim=True
    if not display_defaults.get_value(plot_name,'y min') is None:
        ylow = display_defaults.get_value(plot_name,'y min')
        doSetLim=True
    if doSetLim:
        ax.set_ylim([ylow, yhigh])
    if yunits is not None:
        if  display_defaults.get_value(plot_name,'no_title'):
            ax.set_ylabel(yunits)
        else:
            ax.set_ylabel(ylabel+' (' + yunits + ')')
        ax.yaxis.set_units(yunits)
    else:
        ax.set_ylabel(ylabel)
   

    if times[-1]-times[0] > timedelta(days=3.0) : 
        ax.set_xlabel('Date')
    else:
        ax.set_xlabel('Times (UT)')
    try:
        if times[0].strftime('%d') == times[-1].strftime('%d'):
            time_str = times[0].strftime('%d-%b-%Y')
        else:
            time_str = times[0].strftime('%d-%b')\
              +' --> ' + times[-1].strftime('%d-%b-%Y')
    except:
            try:
             time_str = times.strftime('%d-%b-%Y')
            except:
             time_str = 'UNKNOWN TIME'
    
    if display_defaults.get_value(plot_name,'no_title'):
        ax.set_title(ylabel)
    else:
        ax.set_title(instrument + ' '+title +'  '
                  + time_str)
    if not display_defaults.get_value(plot_name,'no_legend') and legend_list is not None and len(legend_list)>0:
        fontP = FontProperties()
        fontP.set_size('small')
        #legend([plot1], "title", prop = fontP)
        ax.legend(legend_list,loc=legend_position,prop = fontP)
    return figs

def make_xaxis( times, ax):
    ax.set_xlim([times[0], times[-1]])
    [major_locator, minor_locator, major_fmt] = \
        locate_xtime_ticks(times)
    ax.xaxis.set_major_locator(major_locator)
    ax.xaxis.set_major_formatter(major_fmt)
    if times[-1]-times[0] < timedelta(days=1.0) :
        ax.xaxis.set_minor_locator(minor_locator)
    ax.grid(True)
    if times[-1]-times[0] > timedelta(days=3.0) :
        ax.set_xlabel('Date')
    else:    
        ax.set_xlabel('Time (UT)')
    return
 

def graph_setup(f, display_defaults,plot_title,plot_name,clear_fig):

    subplotval=0   
    #f.set_size_inches(display_defaults.get_size('image_size'),forward=True)
    if 'subplot' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
        subplotval=display_defaults.get_value(plot_name,'subplot')
        idx=subplotval%10
        rows=subplotval/100
        cols=(subplotval/10)%10
        minbotrow=(rows-1)*cols
        needXLabel=(idx>minbotrow)
    
    ax=None
    if clear_fig or f.isNewFigure:
        if subplotval>0:
            if not hasattr(f,'subplot_axes'):
                setattr(f,'subplot_axes',dict())
            if subplotval in f.subplot_axes:
                ax=f.subplot_axes[subplotval]
                ax.cla()
                ax=None
        else:
            if hasattr(f,'subplot_axes'):
                for x,a in f.subplot_axes.items():
                    f.delaxes(a)
                setattr(f,'subplot_axes',dict())
            f.clf()
            #f.set_dpi(int(display_defaults.get_value('image_pixels_per_inch','ppi')))
            f.canvas.set_window_title('Fig ' + str(f.number)
                + '       ' +plot_title)
            #f.add_axes([.08,.05,.96,.95])
            ax = f.subplots_adjust(
                top=.9,
                bottom=.1,
                left=.1,
                right=.9,
                hspace=0,
                wspace=0,
                )
    if ax is None:
        if subplotval>0:
            if not hasattr(f,'subplot_axes'):
                setattr(f,'subplot_axes',dict())
            if subplotval in f.subplot_axes:
                ax=f.subplot_axes[subplotval]
                if not ax in f.axes:
                    ax=None
            if ax is None:
                if 'rect' in display_defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
                    rect=display_defaults.get_value(plot_name,'rect')
                    ax=f.add_axes(rect)
                else:
                    ax=f.add_subplot(subplotval)
                f.subplot_axes[subplotval]=ax
            else:
                f.sca(ax)
        else:
            pass #nothing to do
    return

def asSeconds(T,t0=None,filler=np.NAN,returnBasetime=False):
    if t0 is None:
        for x in T:
            if x is not None:
                t0 = x
                break
    ret = []
    for x in T:
        if x is None:
            ret.append(filler)
        else:
            ret.append((x-t0).total_seconds())
    if returnBasetime:
        return ret,t0
    return ret

def defined(x):
    mask=np.zeros_like(x,dtype='bool')
    rmask=mask.ravel()
    for i,v in enumerate(x.ravel()):
        rmask[i] = v is not None
    return mask

def definite(x):
    mask=np.zeros_like(x,dtype='bool')
    rmask=mask.ravel()
    for i,v in enumerate(x.ravel()):
        rmask[i] = np.isfinite(v)
    return mask

def rti_fig(plot_name,instrument,image_array, times, altitudes, title_str, units,
    mask, defaults,figs,fig_name=None,alt_units=None,alt_label=None,homogenize=None):
    """Plots an altitude vs. time image
       images may be plotted with log scaling, linear scaling, or as integer class members
       depending on the parameter defaults.log_linear which may be 'log' | 'linear' | 'cl'
       image array = array to display with dimensions image_array(ntimes,nalts)
                     if defaults.log_linear == 'cl', this must consist of intergers
                     less than len(units)
       times       = array of python datetimes, times(ntimes)
       altitudes   = array of altitudes in meters, altitudes(nalts)
       title_str   = string array for title bar
       cscale	   = min_max values to plot for image_array
                     if defaults.log_linear == 'cl', this is ignored
       units       = units string for color bar
                     
                     if defaults.log_linear == 'cl', this is a list of class names to
                        appear on color bar e.g. with 3 classes units =['none','water','ice']
                        in this mode, any value outside of [0,len-1] is trimmed to the closest extent (0 or len-1)

                     to select colors and values explicitly, use a list of gt.ClassColormapEntry
                        example: [ gt.ClassColormapEntry(0,'none','black'),
                                gt.ClassColormapEntry(1,'water','blue'),
                                gt.ClassColormapEntry(2,'ice','white') ]
                        values can be in any order. classes will be presented on the colorbar in the order provided.
                        If you don't care about the color, set it to None.
                        Any value in the image_array parameter that isn't given in the list default to the value of the first given
       mask        = logical array with same dimensions as image_array,
                     False pixels set to grey
                     if defaults.log_linear == 'cl', this is ignored
       defaults    = json_config instance
       figs         = figure number
       
    """

    import time
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from datetime import datetime
    from matplotlib.font_manager import FontProperties
    import matplotlib.image as imgplot
    import numpy as np
    import matplotlib.image as mpimg

    amask=definite(altitudes)
    tmask=defined(times)
    times=times[tmask]
    altitudes=altitudes[amask]
    tmask=tmask[:,np.newaxis]
    amask=np.transpose(amask[:,np.newaxis])
    newsize=(np.sum(tmask),np.sum(amask))
    #print 'rti shapes',image_array.shape,tmask.shape,amask.shape,(tmask*amask).shape,newsize
    image_array=image_array[tmask*amask].reshape(newsize)
    if mask is not None:
        mask=mask[tmask*amask].reshape(newsize)
    #print 'rti now',image_array.shape

    if alt_units is None:
        alt_units='km'
        altitudes = altitudes/1000.0

    [ntimes, nalts] = image_array.shape

    if ntimes==0 or nalts==0 or times.size==0:
        return

    dtcheck=np.diff(np.array([(x-times[0]).total_seconds() for x in times]))
    avgdt=np.median(dtcheck)
    ddt=np.abs(dtcheck-avgdt)
    dacheck=np.diff(altitudes)
    avgda=np.median(dacheck)
    dda=np.abs(dacheck-avgda)
    reqdif=0.01
    if np.any(ddt>(reqdif*avgdt)) or np.any(dda>(reqdif*avgda)):
        message='rti_fig requires regular time and altitude axes. Threshhold of %.1f%% from da %f %s and dt %f sec Check source!' % (reqdif*100.0,avgda,alt_units,avgdt)
        if homogenize is None:
            homogenize=True
            print 'WARNING: '+message
        if not homogenize:
            print 'altitudes =',altitudes
            print 'da =',dacheck
            print 'dda =',dda
            print 'times =',times
            print 'dt =',dtcheck
            print 'ddt =',ddt
            print 'problem in time =',np.where(ddt>(reqdif*avgdt))
            print 'problem in alt =',np.where(dda>(reqdif*avgda))
            raise RuntimeError(message)
        print 'INTERPOLATING 2D Image for plot '+plot_name
        timesrc=np.array([(x-times[0]).total_seconds() for x in times])
        timedest=np.array([(avgdt*x) for x in range(int(timesrc[-1]/avgdt)+1)])
        newtimes=np.array([(times[0]+timedelta(seconds=x)) for x in timedest])
        #newtimes=type(times)(newtimes)
        newalts=np.array([(altitudes[0]+x*avgda) for x in range(int(((altitudes[-1]-altitudes[0]))/avgda)+1)])
        #newalts=type(altitudes)(newalts)
        arr=image_array
        if False:
            from scipy import interpolate
            func=interpolate.interp2d( altitudes, timesrc, arr, kind='linear',fill_value=None,copy=False)
        else:
            from lg_base.core.array_utils import nearest_nd
            func=nearest_nd([timesrc,altitudes],arr,maxddims=[2*avgdt,2*avgda])
        #print 'interpolating dtimes',timesrc,'to',timedest
        #print 'interpolating alts',altitudes,'to',newalts
        newimage_array=func(timedest,newalts,assume_sorted=True)
        #newimage_array=type(image_array)(newimage_array)
        if mask is not None:
            if False:
                from scipy import interpolate
                func=interpolate.interp2d( altitudes, timesrc, mask, kind='linear',fill_value=None,copy=False)
            else:
                from lg_base.core.array_utils import nearest_nd
                func=nearest_nd([timesrc,altitudes],mask,maxddims=[2*avgdt,2*avgda])
            newmask=func(timedest,newalts,assume_sorted=True)
            mask=newmask

        image_array=newimage_array
        times=newtimes
        altitudes=newalts
        [ntimes, nalts] = image_array.shape

    start_time_str = times[0].strftime('%d-%b-%Y')

    plt.rc('font', size=int(defaults.get_value('image_font_size',"font_points",require=True)))
    
    scaletype =defaults.get_value(plot_name,"log_linear_class",missing=None)
    #check for obsolete name--values for plot type are now: 'log'|'linear'|'class'
    #also accepts: 'log'|'linear'|'cl'
    if scaletype == None:
        scaletype = defaults.get_value(plot_name,"log_linear")
    if scaletype == 'cl':
        scaletype = 'class'
    subplotval=0
    needXLabel=True
    if fig_name is not None:
        pass
    elif 'figure' in defaults.get_labels(plot_name):
        fig_name=defaults.get_value(plot_name,'figure')#to override where this ends up being plotted
    else:
        fig_name=plot_name
    f=figs.figure(fig_name,dpi=int(defaults.get_value('image_pixels_per_inch','ppi',require=True)),figsize=defaults.get_size('image_size'))
    if 'subplot' in defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
        subplotval=defaults.get_value(plot_name,'subplot')
        idx=subplotval%10
        rows=subplotval/100
        cols=(subplotval/10)%10
        minbotrow=(rows-1)*cols
        needXLabel=(idx>minbotrow)
    image_size = defaults.get_size('image_size')
    #f.set_size_inches(image_size, forward=True)
    #f.set_dpi(int(defaults.get_value('image_pixels_per_inch','ppi')))
    
    if subplotval>0:
        if 'rect' in defaults.get_labels(plot_name):#for multiple plots on one figure. all but the first use 'figure' too
            rect=defaults.get_value(plot_name,'rect')
            f.add_axes(rect)
        else:
            f.add_subplot(subplotval)
    else:
        f.add_axes([.08,.05,.96,.95])
        f.canvas.set_window_title('Fig ' + str(f.number)
                + '        '+title_str)
    cscale=np.zeros(2)
    if scaletype in ('log','linear'):
        try:
            cscale[0] = float(defaults.get_value(plot_name,'lo_color_lmt',require=True))
        except KeyError:
            print 'WARNING: missing lo_color_lmt from figure '+plot_name
            cscale[0] = np.nanmin(image_array)
        try:
            cscale[1] = float(defaults.get_value(plot_name,'hi_color_lmt',require=True))
        except KeyError:
            print 'WARNING: missing hi_color_lmt from figure '+plot_name
            cscale[1] = np.nanmax(image_array)
        img=image_array.copy()
        img[~np.isfinite(img)]=np.NaN

        np.seterr(invalid = 'ignore')
        img[img < cscale[0]] = cscale[0]
        img[img > cscale[1]] = cscale[1]
        np.seterr(invalid = 'warn')
        
        if scaletype=='log':
            img = np.log10(img)
            cscale = np.log10([cscale[0],cscale[1]])
        if plot_name == 'gray_masked_image':
            print 'modifing color scale for gray masked image'
            units = defaults.get_value('gray_masked_image','units')
            #move unmasked data to lower half of color scale
            img[mask==1] = (img[mask==1]-cscale[0])/2.0 +cscale[0]
            #move masked data to upper half of color scale
            img[mask==0] = (img[mask==0]-cscale[0])/2.0 + (cscale[0]+cscale[1])/2.0
    elif scaletype == 'class': # integer class plot
        cscale[0] = -0.5
        cscale[1] = len(units)-.5
        img=np.zeros(image_array.shape,dtype='int16')
        cllabels=[]
        clcolors=[]
        defaultcolors=['white','red','orange','yellow','lime','blue','magenta','violet','black']
        if isinstance(units[0],ClassColormapEntry):
            for i,classentry in enumerate(units):
                img[image_array==classentry.value]=i
                cllabels.append(classentry.name)
                if classentry.color is None:
                    clcolors.append(defaultcolors[i%len(defaultcolors)])
                else:
                    clcolors.append(classentry.color)
            units=cllabels
        else:
            for i in range(len(units)):
                clcolors.append(defaultcolors[i%len(defaultcolors)])
    else:
        print
        print
        print 'failed to provide "log_linear_class" value for plot------'+ plot_name
        print
        print
        return
    if defaults.enabled('mask_image') == 1 and not mask is None \
           and not (plot_name == 'gray_masked_image'):
     
    #if defaults.enabled('mask_image') == 1 and not mask is None \
    #       and not (scaletype =='class' or plot_name == 'gray_masked_image'):
        mask = np.bitwise_and(mask.astype('uint16'),np.ones_like(mask.astype('uint16')))
        try:
            img[mask == 0] = np.NaN
            img[np.isnan(image_array)] = np.NaN
        except:
            img[mask == 0] = 0
           
    # expand vertical to fit canvas for small number of vertical pixels

    real_extents=(mdate.date2num(times[0]),mdate.date2num(times[-1]),altitudes[0],altitudes[-1])
    pixel_aspect_ratio = np.float( real_extents[3]-real_extents[2]) \
        / np.float(real_extents[1]-real_extents[0])
    if pixel_aspect_ratio == 0: 
        pixel_aspect_ratio = 1
    canvas_aspect_ratio = np.float(image_size[1])/np.float(image_size[0])
   
    aspect_ratio=0.8*canvas_aspect_ratio/pixel_aspect_ratio
   
    ax=f.gca()
    origintype='lower'
    if 'fliprti' in os.getenv("GRAPHICSHACKS","").split(','):
        origintype=None
    if 'forcefulfliprti' in os.getenv("GRAPHICSHACKS","").split(','):
        img=np.fliplr(img)
    # plot image
    imgplot = ax.imshow(np.transpose(img), aspect=aspect_ratio,
                         origin=origintype,extent=real_extents)

    # imgplot=f.gca().imshow(np.transpose(img),aspect=1,origin='upper')
    # set color limits

    imgplot.set_clim(cscale)
    ax.set_title(title_str + '   ' + start_time_str)
    ax.yaxis.set_units(alt_units)
    if alt_label is not None:
        if len(alt_label)>0:
            ax.set_ylabel(alt_label)
    elif defaults.enabled('image_altitude_above_ground_level'):
        ax.set_ylabel('Altitude--agl ('+alt_units+')')
    else:
        ax.set_ylabel('Altitude--msl ('+alt_units+')')
    
    make_xaxis(times,ax)

    # modify color map to show missing data as black and masked as grey


    myjetcm = cm.jet
    myjetcm.set_bad('k', alpha=1)
    myjetcm.set_under('k',alpha=1)
    imgplot.set_cmap(myjetcm)


    #if gray scale masking
    #split color scale in half--lower half color, upper half gray scale
    if plot_name == 'gray_masked_image' and (scaletype == 'linear' or scaletype == 'log'):
        print 'splitting color scale for gray masked plot'
        g_mask_cm = np.zeros((256,4))
        for i in range(256):
           if i==2*(i/2):
              g_mask_cm[i/2,:] = cm.jet(i)[:]
           else: 
              g_mask_cm[i/2+128,:] = cm.gray(i)[:]
        
        #implement the new color scale
        gray_mask_cm = cm.colors.ListedColormap(g_mask_cm,name=u'from_list') 
        imgplot.set_cmap(gray_mask_cm)

       
        
    
    if scaletype == 'class':
       cl_cmap = colors.ListedColormap(clcolors,'cl_colormap')
       imgplot.set_cmap(cl_cmap)
    

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(1)


    if 1:
        f.gca().grid(True)  # ,linestyle='-')
        if needXLabel:
            if times[-1]-times[0] > timedelta(days=3.0) :            
               f.gca().set_xlabel('Date')
            else:
               f.gca().set_xlabel('Time (UT)')

        # By default, use automatic tick location selection and text formatting
        # for our colorbar

        cb_ticks = None
        cb_format = None
        cb_bounds = None

        # For log plots, set tick locations and formatting explicitly
    
        if scaletype == 'log':
            min_tick = np.int(np.ceil(cscale[0]))
            max_tick = np.int(np.floor(cscale[1]))
            cb_ticks = range(min_tick, max_tick + 1)
            cb_format= expLabels()
        elif scaletype == 'class':
            cb_ticks = range(len(units))
            cb_format= listLabels(units)
    if not hasattr(f,'fig_colorbars') or subplotval not in f.fig_colorbars :
        if not hasattr(f,'fig_colorbars'):
            setattr(f,'fig_colorbars',{})
        cb = f.colorbar(imgplot, ax=f.gca(),shrink=0.7,aspect=15, ticks=cb_ticks,
                          format=cb_format, boundaries=cb_bounds)
        f.fig_colorbars[subplotval]=cb
        if not scaletype == 'class':
            cb.set_label(units)
        #for graph of classes such as mask type
        #if scaletype == 'cl':
        #    cb.ax.yaxis.set_major_formatter(listLabels(units))
                #['no mask','lock','m_lost','cloud','m_lost,i2','cloud,m_lost','all'])
        #if scaletype == 'log':
        #    cb.ax.yaxis.set_major_formatter(expLabels())        
    else:
        print 'COLORBAR for image '+plot_name+' is stale. might not be correct'
    #f.set_size_inches(f.get_size_inches(),forward=True)
    
    # print dir(plt.cm)
    # print dir(plt.cm.get_cmap(plt.cm.Blues,255))
    # print colors
    # print 'cmap',plt.cm.ScalarMappable.get_cmap(plt.jet)


def locate_xtime_ticks(times):
    start = times[0]
    end = times[-1]

    # N==number of ticklabels per hour
    diffseconds=(end-start).total_seconds()
    diffminutes=diffseconds/60.0
    diffhours=diffminutes/60.0
    diffdays=(end-start).days

    if diffdays >= 7:  # plot longer than 7 days
        major_locator = mdate.DayLocator()
        #major_locator = HourLocator(range(0, 25, 12))
        minor_locator = mdate.HourLocator(range(0, 25, 1))
        major_fmt = mdate.DateFormatter('%d')
    elif diffdays >= 3:  # plot longer than 3 days
        major_locator = mdate.DayLocator()
        #major_locator = HourLocator(range(0, 25, 12))
        minor_locator = mdate.HourLocator(range(0, 25, 1))
        major_fmt = mdate.DateFormatter('%m/%d')
    elif diffdays >= 1:  # plot longer than 1 day
        major_locator = mdate.HourLocator(range(0, 25, 4))
        minor_locator = mdate.HourLocator(range(0, 25, 1))
        major_fmt = mdate.DateFormatter('%H')
    elif diffhours >= 12.0:

                                 # plot longer than 12 hours

        major_locator = mdate.HourLocator(range(0, 25, 2))
        minor_locator = mdate.MinuteLocator(range(0, 61, 15))
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffhours >= 4.0:

                                 # plot longer than 4 hours

        major_locator = mdate.HourLocator(range(0, 25, 1))
        minor_locator = mdate.MinuteLocator(range(0, 61, 30))
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffhours >= 2.0:

                                 # plot longer than 2 hours
       # 30 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 30))
        minor_locator = mdate.MinuteLocator(range(0, 61, 5))
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffhours >= 1.0:

                                 # plot longer than 1 hour
       # 15 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 15))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffminutes >= 30.0:

                                 # plot longer than 1/2 hour
       # 10 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 10))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffminutes >= 15.0:

                                 # plot longer than 15 min
       # 5 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 5))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffminutes >= 10.0:

                                 # plot longer than 10 min
       # 2 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 2))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffminutes >= 5.0:

    # elif (end-start)>=1/288.0:   # plot longer than 5 min
                                 # plot longer than  mi
       # 1 min tick labels

        major_locator = mdate.MinuteLocator(range(0, 61, 1))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M')
    elif diffminutes >= 2.0:

                                 # plot longer than 2 min
       # 30 sec tick labels

        major_locator = mdate.SecondLocator(range(0, 61, 30))
        minor_locator = mdate.SecondLocator()
        major_fmt = mdate.DateFormatter('%H:%M:%S')
    elif diffminutes >= 1.0:
        
        # plot longer than 1 min
        # 15 sec tick lables
        
        major_locator = mdate.SecondLocator(range(0, 61, 15))
        minor_locator = mdate.SecondLocator()
        major_fmt = mdate.DateFormatter('%H:%M:%S')

    elif diffseconds >= 30.0:
        
        # plot longer than 30 sec
        # 10 sec tick labels

        major_locator = mdate.SecondLocator(range(0, 61, 10))
        minor_locator = mdate.SecondLocator()
        major_fmt = mdate.DateFormatter('%H:%M:%S')
    else:

                                 # 5 sec tick labels

        major_locator = mdate.SecondLocator(range(0, 61, 5))
        minor_locator = mdate.MinuteLocator()
        major_fmt = mdate.DateFormatter('%H:%M:%S')
    
    
    return (major_locator, minor_locator, major_fmt)

def plot_2d_histogram(plot_name,instrument,times,x_array,y_array,line_type,legend_list,legend_position,xlabel_str\
                     ,x_units,ylabel_str,y_units,title_str,defaults,figs):
    """plot 2-dimensional histogram, first variable ploted as 2d histogram
       optional additonal variables as line plots overploted on histogram.
           plot_2d_histogram(plot_name                #'string'
                          ,instrument                 #'instrument name'
                          ,rs.rs_mean.times              # 
                          ,[x variable array, optional added x  arrays]
                          ,[y variable array, optional added y  arrays]
                          ,['line type']              #line type added vars, ['.r','-g']
                          ,[legend_list]              #legend for optional added plot lines
                          ,legend_position            #legend position, eg. 'upper left'
                          ,xlabel                     #xaxis label
                          ,x_units                    #unit on xaxis
                          ,ylabel                     #yaxis label
                          ,y_units                    #units on yaxis
                          ,title                      #plot title string
                          ,display_defaults           #display defaults .json file
                          ,figs) """
    
    if 1:
        f=figs.figure(plot_name,figsize=defaults.get_size('image_size'))
        f.canvas.set_window_title('Fig ' + str(f.number)
                + '        '+title_str)
        #f.set_size_inches(defaults.get_size('image_size'),forward=True)
   
        number_ybins=float(defaults.get_value(plot_name,'number_ybins',require=True))
        number_xbins=float(defaults.get_value(plot_name,'number_xbins',require=True))
        x_log_linear=defaults.get_value(plot_name,'x_log_linear')
        y_log_linear=defaults.get_value(plot_name,'y_log_linear')
        y_min =defaults.get_value(plot_name,'y min')
        y_max =defaults.get_value(plot_name,'y max')
        x_min =defaults.get_value(plot_name,'x min')
        x_max =defaults.get_value(plot_name,'x max')
        color_log_linear = defaults.get_value(plot_name,'color_log_linear')

        #copy arrays so that changes don't effect other routines using the arrays
        x_ar=x_array[0].copy()
        y_ar=y_array[0].copy()

        #flatten the arrays
        y_ar=np.ravel(y_ar)
        x_ar=np.ravel(x_ar)
        
        #move NaN's outside range of histogram
        #this will place them in bin[0,0] of the histogram
        y_ar[np.isnan(y_ar)]=1e-36
        x_ar[np.isnan(x_ar)]=1e-36
       
        if x_log_linear == 'log':
            assert(x_min>0)
            x_min= np.log10(x_min)
            x_max= np.log10(x_max)
            x_edges=np.logspace(x_min,x_max,number_xbins)
        else:
            x_edges=np.arange(x_min,x_max,(x_max-x_min)/number_xbins)

        if y_log_linear == 'log':
            assert(y_min>0)
            y_min= np.log10(y_min)
            y_max= np.log10(y_max)
            y_edges=np.logspace(y_min,y_max,number_ybins)
        else:
            y_edges=np.arange(y_min,y_max,(y_max-y_min)/number_ybins)
       
        hist_2d,xedges,yedges=np.histogram2d(y_ar,x_ar \
                ,bins=(y_edges,x_edges),normed=False)
        
        #remove the 0,0 element of the histogram because it will contain all
        #of the NaN's
        hist_2d[0,0]=0
                        
        if color_log_linear == 'log':
            hist_2d[hist_2d <= 0] = np.NaN
            histplot = f.gca().imshow(np.log10(hist_2d),origin='lower',aspect='auto'
                ,extent=(x_min,x_max,y_min,y_max))
                #,norm=colors.Normalize(vmin = 0.0,vmax=np.max(np.log10(hist_2d[1:,1:]))))
        else:
            histplot = f.gca().imshow(hist_2d,origin='lower',aspect='auto'
               ,extent=(x_min,x_max,y_min,y_max))
               #,norm=colors.Normalize(vmin = 0.0,vmax=np.max(np.log10(hist_2d[1:,1:]))))
        
        #add over plots, these are stored in x[i],y[i] starting with i=1, x[0],y[0] has histogram data 
        if len(x_array) > 1:
            for i in range(len(x_array)-1):
                if x_log_linear == 'log' and y_log_linear == 'log':
                    f.gca().plot(np.log10(x_array[i+1]),np.log10(y_array[i+1]) ,line_type[i])
                elif x_log_linear == 'log':          
                    f.gca().plot(np.log10(x_array[i+1]),y_array[i+1] ,line_type[i])
                elif  y_log_linear == 'log':       
                    f.gca().plot(x_array[i+1],np.log10(y_array[i+1]) ,line_type[i])
                else:
                    f.gca().plot(x_array[i+1],y_array[i+1] ,line_type[i])
                
        ax = f.gca()
       
        ax.grid(True)
       
        ax.set_xlim([x_min,x_max])
        ax.set_ylim([y_min,y_max])
        

        if x_log_linear == 'log':    
            ax.xaxis.set_major_formatter(logLabels('%6.1e'))
            
        if y_log_linear == 'log':       
            ax.yaxis.set_major_formatter(logLabels('%6.1e'))        
            
        ax.yaxis.grid(color='gray',linestyle='-')
        ax.xaxis.grid(color='gray',linestyle='-')
        time_str=times[0].strftime('%d-%b-%Y %H:%M') \
                 + '--->'+ times[-1].strftime('%H:%M')
        ax.set_title(instrument+' '+title_str+' '+time_str)     
        ax.set_ylabel(ylabel_str+' ('+y_units+')')
        ax.yaxis.set_units(y_units)
        ax.set_xlabel(xlabel_str+' ('+x_units+')')
        ax.xaxis.set_units(x_units)
        if not legend_list == None and len(legend_list) > 0:
            ax.legend(legend_list,loc=legend_position)

        
