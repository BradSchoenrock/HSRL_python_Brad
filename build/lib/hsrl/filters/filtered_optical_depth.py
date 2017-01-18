import numpy as np
import lg_base.core.decoratortools as nt

def filtered_optical_depth(molecular,beta_a,window,order):


    if not (window%2) == 0:
        window = window +1
        print 'savitzky_golay: requires odd window length--adding one'

    data_len = len(y)
    half_win = (window-1)/2
    
    #reflect data ends for padding
    yy=zeros(len(y)+window)
    yy[0:half_win] = y[half_win:0:-1]
    yy[range(data_len-half_win,data_len)] = range(data_len,(data_len-half_win),-1)
    xx = range(data_len)
    
    for i in range(data_len):

        ylocal = yy[i-half_win:i+half-win]  
        pc = np.polyfit(xx,yy,order)
        center_slope = np.polyval(pc[0:order],half_win+1)
        print i, pc, center_slope
        if center_slope < 0:
            pass #THIS CODE WAS OBVIOUSLY NEVER RUN FIXME OR DELETE

    return filtered_od        
