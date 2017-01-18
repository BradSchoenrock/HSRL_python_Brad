'''
International Standard Atmosphere utilities
'''

import math

# Acceleration due to gravity at Earth's surface
_LITTLE_G_ = 9.80655      # m/s^2

# Gas constant for air
_R_ = 287.053             # J/kg K

# The standard atmosphere has 7 layers. Create an array of 7 empty dictionaries
__layers__ = [{}, {}, {}, {}, {}, {}, {}]

# Surface layer starts at the 0 m with temperature 15 C, 
# pressure 1013.25 hPa, and temperature lapse rate of -0.0065 K/m
__layers__[0]['height'] = 0.0                   # m
__layers__[0]['lapse_rate'] = -0.0065           # K/m
__layers__[0]['t'] = 288.15                     # K (= 15 C)
__layers__[0]['p'] = 1013.25                    # hPa

# Layer starting at 11000 m has lapse rate of 0.0 K/m
__layers__[1]['height'] = 11000.0
__layers__[1]['lapse_rate'] = 0.0

# Layer starting at 20000 m has lapse rate of 0.001 K/m
__layers__[2]['height'] = 20000.0
__layers__[2]['lapse_rate'] = 0.001

# Layer starting at 32000 m has lapse rate of 0.0028 K/m
__layers__[3]['height'] = 32000.0
__layers__[3]['lapse_rate'] = 0.0028

# Layer starting at 47000 m has lapse rate of 0.0 K/m
__layers__[4]['height'] = 47000.0
__layers__[4]['lapse_rate'] = 0.0

# Layer starting at 51000 m has lapse rate of -0.0028 K/m 
__layers__[5]['height'] = 51000.0
__layers__[5]['lapse_rate'] = -0.0028

# Layer starting at 71000 m has lapse rate of -0.002 K/m
__layers__[6]['height'] = 71000.0
__layers__[6]['lapse_rate'] = -0.002

# The top of the standard atmosphere is 84852 m
StdAtmosTop = 84852.0

# Calculate the temperature and pressure at the base of each level above the
# surface
for i in range(1, 7):
    hdiff = __layers__[i]['height'] - __layers__[i-1]['height']
    # temperature at the base of this layer
    __layers__[i]['t'] = __layers__[i-1]['t'] + \
        __layers__[i-1]['lapse_rate'] * (hdiff)
    # pressure at the base of this layer
    __layers__[i]['p'] = __layers__[i-1]['p']
    if (__layers__[i-1]['lapse_rate'] == 0.0):
        __layers__[i]['p'] *= math.exp((-_LITTLE_G_ / _R_) * hdiff / __layers__[i-1]['t'])
    else:
        __layers__[i]['p'] *= math.exp(math.log(__layers__[i-1]['t'] / __layers__[i]['t']) * 
                                   (_LITTLE_G_ / _R_) / __layers__[i-1]['lapse_rate'])
        
def LayerBoundHeights():
    '''
    Return an array containing the boundary heights between standard atmosphere 
    layers, in m.
    
    The first and last values of the array define the bottom and top of the
    standard atmosphere, respectively.
    '''
    heights = []
    for layer in __layers__:
        heights.append(layer['height'])
        
    heights.append(StdAtmosTop)
    return heights

def LayerBoundPressures():
    '''
    Return an array containing the boundary pressures between standard 
    atmosphere layers, in hPa.
    
    The first and last values of the array give the bottom and top level
    pressures of the standard atmosphere, respectively.
    '''
    pressures = []
    for h in LayerBoundHeights():
        pressures.append(HeightToPres(h))
        
    return pressures

def LayerBoundTemps():
    '''
    Return an array containing the boundary temperatures between standard 
    atmosphere layers, in K.
    
    The first and last values of the array give the bottom and top level
    temperatures of the standard atmosphere, respectively.
    '''
    temps = []
    for h in LayerBoundHeights():
        temps.append(HeightToTemp(h))
        
    return temps

def PresToHeight(pres):
    '''
    Return standard atmosphere height in m for a given pressure in hPa.
    
    Return standard atmosphere height in m for a given pressure in hPa.
    NaN is returned if the pressure is outside the limits of the standard
    atmosphere.
    '''
    # Undefined at pressures above 1013.25 hPa
    if (pres > __layers__[0]['p']):
        return float('nan')
    
    # Find the layer containing this pressure
    for i in range(6, -1, -1):
        if (pres <= __layers__[i]['p']):
            break
    # Find the height for the given pressure
    h0 = __layers__[i]['height']
    p0 = __layers__[i]['p']
    t0 = __layers__[i]['t']
    lapse_rate = __layers__[i]['lapse_rate']
    
    h = h0
    if (lapse_rate == 0.0):
        h += -math.log(pres / p0) * t0 / (_LITTLE_G_ / _R_)
    else:
        h += (__layers__[i]['t'] / 
              math.pow((pres / p0), (lapse_rate / (_LITTLE_G_ / _R_))) - t0) / lapse_rate             
             
    # Return the calculated height, or NaN if the height is above the top of
    # the standard atmosphere
    if (h > StdAtmosTop):
        return float('nan')
    else:
        return h

#
# Return standard atmosphere temperature in K for a given height in m. NaN is 
# returned if the given height is outside the bounds of the standard atmosphere.
#
def HeightToTemp(h):
    '''
    Return standard atmosphere temperature in K for a given height in m.
    
    Return standard atmosphere temperature in K for a given height in m.
    NaN is returned if the height is outside the limits of the standard
    atmosphere.
    '''
    # Standard atmosphere is undefined below 0 m and above StdAtmosTop
    if ((h < 0.0) or (h > StdAtmosTop)):
        return float('nan')
    
    # Find the layer containing this height
    for i in range(6, -1, -1):
        if (h > __layers__[i]['height']):
            break
        
    # Calculate temperature at this height
    hdiff = h - __layers__[i]['height']
    t = __layers__[i]['t'] + __layers__[i]['lapse_rate'] * hdiff;
    return t
        

def PresToTemp(h):
    '''
    Return standard atmosphere temperature in K for a given pressure in hPa.
    
    Return standard atmosphere temperature in K for a given pressure in hPa.
    NaN is returned if the pressure is outside the limits of the standard
    atmosphere.
    '''
    return HeightToTemp(PresToHeight(h))

def HeightToPres(h):
    '''
    Return standard atmosphere pressure in hPa for a given height in m.
    
    Return standard atmosphere pressure in hPa for a given height in m.
    NaN is returned if the height is outside the limits of the standard
    atmosphere.
    '''
    # Standard atmosphere is undefined below 0 m and above StdAtmosTop
    if ((h < 0.0) or (h > StdAtmosTop)):
        return float('nan')
    
    # Find the layer containing this height
    for i in range(6, -1, -1):
        if (h > __layers__[i]['height']):
            break
    
    # Calculate the associated pressure
    t0 = __layers__[i]['t']  # temperature at bottom of this layer
    p0 = __layers__[i]['p']  # pressure at bottom of this layer
    h0 = __layers__[i]['height']
    lapse_rate = __layers__[i]['lapse_rate']    # lapse rate for this layer

    pres = p0
    hdiff = h - h0
    if (lapse_rate == 0.0):
        pres *= math.exp(-(_LITTLE_G_ / _R_) * hdiff / t0)
    else:
        pres *= math.pow((t0 / (t0 + lapse_rate * hdiff)),
                         (_LITTLE_G_ / _R_) / lapse_rate)
  
    return pres;
    
