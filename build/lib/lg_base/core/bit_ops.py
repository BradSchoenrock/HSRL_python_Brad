import numpy as np
def set_bit(array,trigger_array,bit_number):
    """set_bit(array,trigger_array,bit_number)
       set bit_number in array if trigger == True
       clear bit_number in array if trigger == false
       at the array elements corresponding to elements of
       the trigger array.
       trigger_array must be the same size as the array"""

    shift_array = bit_number*np.ones_like(array)
    clr_mask  = np.left_shift(np.ones_like(array),shift_array)
    clr_mask  = np.invert(clr_mask)
    set_mask = np.left_shift(trigger_array,shift_array)
    #clear the bit position
    array = np.bitwise_and(array,clr_mask)
    #set the bit if corresponding element of trigger array == True
    array = np.bitwise_or(array,set_mask)

    return array

def bit_zero_and(array,trigger_array):
    """and the trigger array which is boolain
       with bit[0] of array"""

    #set all bits except zero bit
    mask = np.invert(np.ones_like(array))
    #set zero bit in mask if trigger_array == True
    mask = np.bitwise_or(mask,trigger_array)

    #and will clear bit zero if either array[0] or trigger
    #is zero
    
    return np.bitwise_and(array,mask)
    
    
    
    
       
    
