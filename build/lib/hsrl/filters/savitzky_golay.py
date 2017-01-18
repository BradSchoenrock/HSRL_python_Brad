#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import scipy, scipy.signal

import lg_base.core.decoratortools as nt

try:
    import parakeet
except ImportError:
    import lg_base.core.dumba as parakeet

def savitzky_golay( y, window_size, order, deriv = 0 ,m=None,retval=False,prepadded=False):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
The Savitzky-Golay filter removes high frequency noise from data.
It has the advantage of preserving the original shape and
features of the signal better than other types of filtering
approaches, such as moving averages techhniques.
Parameters
----------
y : array_like, shape (N,)
the values of the time history of the signal.
window_size : int
the length of the window. Must be an odd integer number.
order : int
the order of the polynomial used in the filtering.
Must be less then `window_size` - 1.
deriv: int
the order of the derivative to compute (default = 0 means only smoothing)
Returns
-------
ys : ndarray, shape (N)
the smoothed signal (or it's n-th derivative).
Notes
-----
The Savitzky-Golay is a type of low-pass filter, particularly
suited for smoothing noisy data. The main idea behind this
approach is to make for each point a least-square fit with a
polynomial of high order over a odd-sized window centered at
the point.
Examples
--------
t = np.linspace(-4, 4, 500)
y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
ysg = savitzky_golay(y, window_size=31, order=4)
import matplotlib.pyplot as plt
plt.plot(t, y, label='Noisy signal')
plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
plt.plot(t, ysg, 'r', label='Filtered signal')
plt.legend()
plt.show()
References
----------
.. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
Data by Simplified Least Squares Procedures. Analytical
Chemistry, 1964, 36 (8), pp 1627-1639.
.. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
Cambridge University Press ISBN-13: 9780521880688
"""
    try:
        window_size = np.abs( np.int( window_size ) )
        order = np.abs( np.int( order ) )
    except ValueError, msg: 
        raise ValueError( "window_size and order have to be of type int" )
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError( "window_size size must be a positive odd number" )
    #if window_size < order + 2:
    if window_size < order + 1:    
        print 'window_size= ',window_size, 'order= ',order
        raise TypeError( "window_size is too small for the polynomials order" )
    if window_size > len(y):
        print 'window_size =', window_size, 'data_length= ' ,  len(y)
        raise ValueError, 'ERROR---Savitzky_golay filter, window size longer than data'
    return __savitzky_golay(y, window_size, order, deriv,m,retval,prepadded)

#@nt.Delayed(parakeet.jit,minversion=dict(parakeet='0.22'))
def __savitzky_golay( y, window_size, order, deriv ,_m,retval,prepadded):
    window_size = np.abs( np.int( window_size ) )
    half_window = ( window_size - 1 ) // 2
    if _m is None:
        window_range = np.arange( -half_window, half_window + 1 )
        order = np.abs( np.int( order ) )
        order_range = np.arange( order + 1 )
        tmp=[]
        for k in window_range:
            tmp2=[]
            for i in order_range:
                tmp2.append(k**i)
            tmp.append(tmp2)
        # precompute coefficients
        b = np.mat(tmp)#np.pow(window_range,order_range)#np.mat( [[k ** i for i in order_range] for k in np.range( -half_window, half_window + 1 )] )
        m = np.linalg.pinv( b ).A[deriv]
    else:
        m=_m
    # pad the signal at the extremes with
    # values taken from the signal itself
    if not prepadded:
        firstvals = y[0] - np.abs( y[1:half_window + 1][::-1] - y[0] )
        lastvals = y[-1] + np.abs( y[-half_window - 1:-1][::-1] - y[-1] )
        y = np.concatenate( ( firstvals, y, lastvals ) )
    if not retval:
        return np.convolve( m, y, mode = 'valid' )
    else:
        return np.convolve( m, y, mode = 'valid' ),m

#@nt.Delayed(parakeet.jit,minversion=dict(parakeet='0.22'))
def __savitzky_golay3( y, window_size, order ):
    return __savitzky_golay(y, window_size, order, 0,None,False,False)

def savitzky_golay_piecewise( xvals, data, kernel = 11, order = 4 ):
    return __savitzky_golay_piecewise( xvals, data, kernel , order )

#@nt.Delayed(parakeet.jit,minversion=dict(parakeet='0.22'))
def __savitzky_golay_piecewise( xvals, data, kernel , order ):
    turnpoint = 0
    last = len( xvals )
    if xvals[1] > xvals[0] : #x is increasing?
        for i in range( 1, last ) : #yes
            if xvals[i] < xvals[i - 1] : #search where x starts to fall
                turnpoint = i
                break
    else: #no, x is decreasing
        for i in range( 1, last ) : #search where it starts to rise
            if xvals[i] > xvals[i - 1] :
                turnpoint = i
                break
    if turnpoint == 0 : #no change in direction of x
        return __savitzky_golay3( data, kernel, order )
    else:
        #smooth the first piece
        firstpart = __savitzky_golay3( data[0:turnpoint], kernel, order )
        #recursively smooth the rest
        rest = __savitzky_golay_piecewise( xvals[turnpoint:], data[turnpoint:], kernel, order )
        return np.concatenate( ( firstpart, rest ) )

#this wrapper function exists because 1) default parameters, 2) exceptions from parameter checking. neither are supported in numba-jit'ed functions
def sgolay2d ( z, window_size, order, derivative = None ):
    """
"""
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2 ) / 2.0

    if window_size % 2 == 0:
        raise ValueError( 'window_size must be odd' )

    if window_size ** 2 < n_terms:
        raise ValueError( 'order is too high for the window size' )
    dval={None:0
          ,'col':1
          ,'row':2
          ,'both':3#compiled version can't change return parameters to accomodate this
          }

    return __sgolay2d( z, window_size, order, dval[derivative])

#@nt.Delayed(parakeet.jit,minversion=dict(parakeet='0.22'))
def __sgolay2d ( z, window_size, order, derivative ):
    # number of terms in the polynomial expression
    #n_terms = ( order + 1 ) * ( order + 2 ) / 2.0

    #if window_size % 2 == 0:
    #    raise ValueError( 'window_size must be odd' )

    #if window_size ** 2 < n_terms:
    #    raise ValueError( 'order is too high for the window size' )

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ ( k - n, n ) for k in range( order + 1 ) for n in range( k + 1 ) ]

    # coordinates of points
    ind = np.arange( -half_size, half_size + 1, dtype = np.float64 )
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1] ).reshape( window_size ** 2, )

    # build matrix of system of equation
    A = np.empty( ( window_size ** 2, len( exps ) ) )
    for i in enumerate( exps ):
        A[:, i[0]] = ( dx ** i[1][0] ) * ( dy ** i[1][1] )

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2 * half_size, z.shape[1] + 2 * half_size
    Z = np.zeros( ( new_shape ) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs( np.flipud( z[1:half_size + 1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs( np.flipud( z[-half_size - 1:-1, :] ) - band )
    # left band
    band = np.tile( z[:, 0].reshape( -1, 1 ), [1, half_size] )
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size + 1] ) - band )
    # right band
    band = np.tile( z[:, -1].reshape( -1, 1 ), [1, half_size] )
    Z[half_size:-half_size, -half_size:] = band + np.abs( np.fliplr( z[:, -half_size - 1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs( np.flipud( np.fliplr( z[1:half_size + 1, 1:half_size + 1] ) ) - band )
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs( np.flipud( np.fliplr( z[-half_size - 1:-1, -half_size - 1:-1] ) ) - band )

    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs( np.flipud( Z[half_size + 1:2 * half_size + 1, -half_size:] ) - band )
    # bottom left corner
    band = Z[-half_size:, half_size].reshape( -1, 1 )
    Z[-half_size:, :half_size] = band - np.abs( np.fliplr( Z[-half_size:, half_size + 1:2 * half_size + 1] ) - band )

    alg = np.linalg.pinv( A )
    # solve system and convolve
    if derivative == 0:#None:
        m = alg[0].reshape( ( window_size, -1 ) )
        return scipy.signal.fftconvolve( Z, m, mode = 'valid' )
    elif derivative == 1:#'col':
        c = alg[1].reshape( ( window_size, -1 ) )
        return scipy.signal.fftconvolve( Z, -c, mode = 'valid' )
    elif derivative == 2:#'row':
        r = alg[2].reshape( ( window_size, -1 ) )
        return scipy.signal.fftconvolve( Z, -r, mode = 'valid' )
    elif derivative == 3:#'both':
        c = alg[1].reshape( ( window_size, -1 ) )
        r = alg[2].reshape( ( window_size, -1 ) )
        return scipy.signal.fftconvolve( Z, -r, mode = 'valid' ), scipy.signal.fftconvolve( Z, -c, mode = 'valid' )
