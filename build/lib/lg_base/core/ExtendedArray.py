import numpy as np

class ExtendedArray(np.ndarray):

    def __new__(subtype, data, info=None, dtype=None, copy=False):
        #        print "__new__ received subtype = '%s', %s, info = %s" % (repr(subtype), type(data), info)
        # Make sure we are working with an array, and copy the data if requested
        subarr = np.array(data, dtype=dtype, copy=copy)

        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)

        # Use the specified 'info' parameter if given
        if info is not None:
            subarr.info = info
        # Otherwise, use data's info attribute if it exists
        elif hasattr(data, 'info'):
                subarr.info = data.info

        # Finally, we must return the newly created object:
        return subarr

    def __array_finalize__(self,obj):
        # We use the getattr method to set a default if 'obj' doesn't have the 'info' attribute
        self.info = getattr(obj, 'info', {})
        # We could have checked first whether self.info was already defined:
        #if not hasattr(self, 'info'):
        #    self.info = getattr(obj, 'info', {})

    def __repr__(self):
        desc="""\
array(data=
  %(data)s)"""
        return desc % {'data': str(self), 'info':self.info }
