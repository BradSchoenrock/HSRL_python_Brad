#!/usr/bin/python
import unittest
import numpy as np
from array_utils import au
#import hsrl.utils.hsrl_array_utils
import types

class TestHsrlArrayUtils(unittest.TestCase):
    def setUp(self): 
        self.t1 = au.T_Array(np.arange(10))
        self.t2 = self.t1 * np.arange(10)
        self.z1 = au.Z_Array(np.arange(10))
        self.tz1 = au.TZ_Array([])
        self.tz_group1 = au.Time_Z_Group(self.t1)
        self.tz_group2 = au.Time_Z_Group(self.t1)
        self.tz_group1.foobar = self.t1
        self.tz_group1.zoo = au.T_Array([])
        self.tz_group2.alpha = self.t1
        self.tz_group2.beta = au.T_Array([])
        self.tz_group2.gamma = 3


    def test_check_id(self):
        self.assertTrue(type(self.t1) == au.T_Array)
        self.assertTrue(type(self.t2) == au.T_Array)
        self.assertTrue(self.tz1.check_type(self.tz1))
        self.assertTrue(self.tz1.check_type2())
        self.assertTrue(self.t1.check_type())
        self.assertTrue(self.tz_group1.check_tarray('foobar'))
        self.assertTrue(self.tz_group1.check_tarray('zoo'))

        self.assertTrue(self.tz_group1.check_tarray2(self.tz_group2, 'alpha'))
        self.assertTrue(self.tz_group1.check_tarray2(self.tz_group2, 'beta'))
        self.assertTrue(au.check_tarray3(self.tz_group2, 'alpha'))
        self.assertFalse(au.check_tarray3(self.tz_group2, 'gamma'))
        self.assertFalse(au.check_tarray3(3, 'gamma'))

if __name__ == '__main__':
    unittest.main()
