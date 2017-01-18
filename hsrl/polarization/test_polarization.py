"""
test_polarization.py

"""
import unittest
import numpy as np
import HSRL_PolCorrection as hp
import csv
from scipy.io import loadmat

polCalTerm = {
"tx_elip": -0.036404364,
"cross_rx_diat": 0.88127382,
"cross_rx_rot": 1.6024794,
"cross_rx_elip": 0.046456464,
"mol_rx_diat": 0.99611489,
"mol_rx_rot": 0.08783443,
"mol_rx_elip": -0.0082449222,
"total_rx_diat": 0.99985812,
"total_rx_rot": 3.2317574,
"total_rx_elip": -0.0086052231,
"qwp_phase_0": 1.5350108,
"qwp_phase_2": -1.08365,
"qwp_rot_0": 0.35391488,
"mirror_phase": 3.3447555,
"mirror_diat": -0.010627414,
"mirror_rot": -1.1491565,
"qwp_rot_offset": -0.041781656,
"qwp_phase_1": 0.095061646,
"qwp_rot_1": -2.5950543,
"Empty_Slot": 1,
"window_phase_rot": -0.0055755222,
"cal_total_depol": 0.98312967,
"cal_molecular_depol": 0.99663812,
"tx_rot": -1.5519834,
"window_phase": -0.073085796,
"window_diat": 0.012753727,
"window_diat_rot": 0.00014232693,
"lowgain_rx_diat": 0.98623903,
"lowgain_rx_rot": 3.256742,
"lowgain_rx_elip": -0.019271303,
}

class TestPolarization(unittest.TestCase):
    def setUp(self):
        self.diatuValidation = loadmat('DiatU_Validation.mat')
        self.VWPValidation = loadmat('VWP_Validation.mat')
        self.RevMVal = loadmat('RevMueller_Validation.mat')
#        self.polCor = loadmat("HSRLPolConfig_InternalVars.mat")
        self.polCor = loadmat("HSRLPolConfigLH_Validation.mat")
        self.polCalTerm = polCalTerm
        
    def test_diatu(self):
        epsilon = 1e-10
        # retrieve inputs and expected outputs
        d_input = self.diatuValidation['D_input']
        d_output= self.diatuValidation['DiatUoutput']
        for i in range(d_input.shape[1]):
            out = hp.DiatU(d_input[:,i])            
            maxdiff = np.max(np.abs(out - d_output[:,:,i]))
            self.assertTrue(maxdiff < epsilon)
        
    def test_VWP(self):
        epsilon = 1e-8
        # retrieve inputs and expected outputs
        theta = self.VWPValidation['theta_input']
        gamma = self.VWPValidation['Gamma_input']
        expected = self.VWPValidation['VWPoutput']
        for i in range(theta.shape[1]):
            out = hp.VWP(theta[0,i], gamma[0,i])
            maxdiff = np.max(np.abs(out - expected[:,:,i]))
            self.assertTrue(maxdiff < epsilon)    
            
    def test_RevMueller(self):
        epsilon = 1e-8
        # retrieve inputs and expected outputs
        gamma = self.RevMVal['RevMueller_input']
        expected = self.RevMVal['RevMueller_output']
        for i in range(gamma.shape[1]):
            out = hp.RevMueller(gamma[:,:,i])
            maxdiff = np.max(np.abs(out - expected[:,:,i]))
            self.assertTrue(maxdiff < epsilon)         

    def internal_test_polarization(self ):
        epsilon = 1e-6
        """test the HSRL polarization routines by comparing with known outputs
        """
#        (Si,DxV, DtV, DmV)= hp.HSRLPolConfig(self.polCalTerm, 
#                                 self.polCor['QWPAngle'][0])
        angles = np.zeros((1))
        dict = hp.HSRLPolConfig(self.polCalTerm, angles)
        for k in dict.keys():
            print 'Comparing ', k
            expected = self.polCor[k]
            out = dict[k]
            maxdiff = np.max(np.abs(out - expected))
            self.assertTrue(maxdiff < epsilon)         
            
        print 'test 1 complete'
        
    def check_array(self, value, expected, epsilon):
        maxdiff = np.max(np.abs(value - expected))
        self.assertTrue(maxdiff < epsilon)         
        
    def test_polarization(self ):
        epsilon = 1e-6
        """test the HSRL polarization routines by comparing with known outputs
        """
        (Si,DxV, DtV, DmV)= hp.HSRLPolConfig(self.polCalTerm, 
                                 self.polCor['thetaIn'][0])
        self.check_array(Si, self.polCor['So'], epsilon) 
        self.check_array(DxV, self.polCor['Dx'], epsilon) 
        self.check_array(DtV, self.polCor['Dt'], epsilon) 
        self.check_array(DmV, self.polCor['Dm'], epsilon)
        print 'test polarization complete'    

if __name__ == '__main__':
    unittest.main()
    
