import numpy as np
from numpy import transpose,trapz
from math import cos, sin
from numpy.linalg import pinv
import lg_base.core.array_utils as hau
from datetime import timedelta
# -*- coding: ascii -*-

class HSRL_PolCorrection:

    def __init__(self, rs_mean, rs_cal, requested_times):
        """
        rs_mean - stores original data
        rs_cal - calibration information
        'requested_times' are bin edges of requested intervals. (40+ second)
        """
        self.rs_mean = rs_mean
        self.rs_cal = rs_cal
        self.requested_times = requested_times

    def computeCorrection(self):
        """
        compute corrected values given a set of rotation angles for the quarter wave plate
        """

        # convert back to radians (squaring conversion factor is temporary fix.
        # somewhere in processing code qwp rotation is being converted to 
        # degrees twice.)
        qw_rot_angle = self.rs_mean.qw_rotation_angle * (np.pi/ 180.0)  
        #qw_rot_angle = self.rs_mean.qw_rotation_angle

        # estimate the previous ('-1'st) angle
        qw_rot_angle = np.append(2* qw_rot_angle[0] - qw_rot_angle[1], qw_rot_angle)

        qw_integ_res = 1 * np.pi / 180.0  #  make this a processing parameter?
        num_angles = len(qw_rot_angle)-1  # Matt: added -1
        measX10 = np.zeros( (num_angles, 10))
        measT10 = np.zeros( (num_angles, 10))
        measX2 = np.zeros( (num_angles, 2))
        measT2 = np.zeros( (num_angles, 2))
        measM2 = np.zeros( (num_angles, 2))

        #debugVar = np.zeros((num_angles,))
        for index in range(len(qw_rot_angle)-1):
            subQWP = np.arange(qw_rot_angle[index], qw_rot_angle[index+1], qw_integ_res)
            #Matt Add - new condition when QWP is not rotating to avoid empty array
            if np.size(subQWP) == 0:
                subQWP = np.array([qw_rot_angle[index]])
            elif np.size(subQWP) == 1:
                subQWP = np.array([subQWP])
                
            [Si,DxV,DtV,DmV]  = HSRLPolConfig(self.rs_cal.pol_cal.data, subQWP)

            # DxV = DxV*0.30
            #  Form 10 element measurement Matrices
            measX10[index,:] = HSRL_FormMeasMatrix(Si,DxV,10)  # 1 x 10
            measT10[index,:] = HSRL_FormMeasMatrix(Si,DtV,10)  # 1 x 10
            #  Form 2 element Measurement Matrices
            measX2[index,:] = HSRL_FormMeasMatrix(Si,DxV,2)  # 1 x 2
            measT2[index,:] = HSRL_FormMeasMatrix(Si,DtV,2) # 1 x 2
            measM2[index,:] = HSRL_FormMeasMatrix(Si,DmV,2) # 1 x 2
            
            #debugVar[index] = num_angles

        # original times are rs_mean.times (Python Datetime)

        if self.requested_times==None:
            n_times=self.rs_mean.times.shape[0]
            bin_edges = self.rs_mean.times.copy()
            if n_times>0:
                bin_edges[n_times]=bin_edges[n_times-1]+timedelta(seconds=.01)
        else:
            n_times=self.requested_times.shape[0]-1
            bin_edges = self.requested_times.copy()
        
        # Matt: Create bin_increments list containing start and stop index of each column
        # Matt: Assumes rs_mean.times contains the supplied data time array (at half second intervals)
        ind = np.arange(self.rs_mean.times.shape[0])

        # make a list of indices to indicate which elements to average for each entry
        # in new times vector
        k = 0
        temp = self.rs_mean.times[ :]
        bin_increments = [];
        while k < bin_edges.shape[0]-1 :
            inxa = ind[temp >= bin_edges[ k]]
            temp2 = self.rs_mean.times[ inxa]
            inx = inxa[temp2 < bin_edges[ k + 1]]
            bin_increments.append([inx[0], inx[-1]])
            k = k + 1
            
        

        cnt = len(bin_increments)
        #TODO - is shape[1] correct???  Matt:  Looks like it is correct.
        nalts = self.rs_mean.cross_pol_counts.shape[1]
        self.counts_inv = np.zeros(  (2*cnt, nalts) )
        self.mol_counts_inv = np.zeros(  (2*cnt, nalts) )
        self.virtual_total_counts = np.zeros( (cnt, nalts))
        self.virtual_cross_counts = np.zeros( (cnt, nalts))
        self.virtual_mol_counts = np.zeros( (cnt, nalts))
        self.f11 = np.zeros( (cnt, nalts))
        self.f12 = np.zeros( (cnt, nalts))
        self.f13 = np.zeros( (cnt, nalts))
        self.f14 = np.zeros( (cnt, nalts))
        self.f22 = np.zeros( (cnt, nalts))
        self.f23 = np.zeros( (cnt, nalts))
        self.f24 = np.zeros( (cnt, nalts))
        self.f33 = np.zeros( (cnt, nalts))
        self.f34 = np.zeros( (cnt, nalts))
        self.f44 = np.zeros( (cnt, nalts))

        for index in range(len(bin_increments)-1):
            # Compute 10 element total scattering matrix  
            Minv_10element = pinv(np.append(
                measX10 [bin_increments[index][0]:bin_increments[index][1],:],
                measT10[bin_increments[index][0]:bin_increments[index][1],:], axis=0) )

            counts_inv = np.append((self.rs_mean.cross_pol_counts[bin_increments[index][0]:bin_increments[index][1],:]), 
                                   (self.rs_mean.combined_counts[bin_increments[index][0]:bin_increments[index][1],:]), axis=0)
            f_10element= np.dot( Minv_10element ,counts_inv)  


            # Compute 2 element total scattering matrix (assumes randomly oriented scatterers)
            Minv_2element = pinv(np.append(measX2 [bin_increments[index][0]:bin_increments[index][1],:],  
                                           measT2[bin_increments[index][0]:bin_increments[index][1],:], axis=0))
            f_2element= np.dot(Minv_2element,counts_inv)

            # Compute Molecular scattering matrix
            mol_counts_inv = (self.rs_mean.molecular_counts[bin_increments[index][0]:bin_increments[index][1],:]) 	
            Minv_mol = pinv(measM2 [bin_increments[index][0]:bin_increments[index][1],:])
            f_mol= np.dot(Minv_mol,mol_counts_inv)

            # Matt:  Transposed indexing on assignment variables
            self.virtual_total_counts[index,:] = (f_2element[1,:])
            self.virtual_cross_counts[index,:] = (f_2element[0,:]-f_2element[1,:])
            #self.virtual_cross_counts[index,0] = debugVar[index]
            self.virtual_mol_counts[index,:] = (f_mol[1,:])
            # save scattering matrix elements
            self.f11[index,:] = f_10element[0,:]
            self.f12[index,:] = f_10element[1,:]
            self.f13[index,:] = f_10element[2,:]
            self.f14[index,:] = f_10element[9,:]
            self.f22[index,:] = f_10element[3,:]
            self.f23[index,:] = f_10element[4,:]
            self.f24[index,:] = f_10element[5,:]
            self.f33[index,:] = f_10element[6,:]
            self.f34[index,:] = f_10element[7,:]
            self.f44[index,:] = f_10element[8,:]

    def applyCorrection(self, rs_mean):
        rs_mean.cross_pol_counts = hau.TZ_Array(self.virtual_cross_counts)
        rs_mean.combined_counts = hau.TZ_Array(self.virtual_total_counts)
        rs_mean.molecular_counts = hau.TZ_Array(self.virtual_mol_counts)
        # insert other computed variables here
        rs_mean.f11 = hau.TZ_Array(self.f11)
        rs_mean.f12 = hau.TZ_Array(self.f12)
        rs_mean.f13 = hau.TZ_Array(self.f13)
        rs_mean.f14 = hau.TZ_Array(self.f14)
        rs_mean.f22 = hau.TZ_Array(self.f22)
        rs_mean.f23 = hau.TZ_Array(self.f23)
        rs_mean.f24 = hau.TZ_Array(self.f24)
        rs_mean.f33 = hau.TZ_Array(self.f33)
        rs_mean.f34 = hau.TZ_Array(self.f34)
        rs_mean.f44 = hau.TZ_Array(self.f44)

# Matlab Code for Measurement Matrix Calculations:

def HSRLPolConfig(x,thQWP):
    """
outputs the polarization configuration vectors of all 4 HSRL channels for a provided
    x - a dictionary containing all polarization calibration terms
    thQWP - a range of QWP angles in radians
returns 
[Si,DxV,DtV,DmV] 

"""


    thQWP = thQWP+x['qwp_rot_offset']  

    DxV = np.zeros((len(thQWP),4))  # N x 4
    DtV = np.zeros((len(thQWP),4))  # N x 4
    DmV = np.zeros((len(thQWP),4))  # N x 4
    Si = np.zeros((4,len(thQWP)))  # 4 x N

    Stx = np.array((1,cos(2*x['tx_elip'])*cos(2*x['tx_rot']),cos(2*x['tx_elip'])*sin(2*x['tx_rot']),sin(2*x['tx_elip'])))  # 4 x 1
    Dx = np.array((1,x['cross_rx_diat']*cos(2*x['cross_rx_rot'])*cos(2*x['cross_rx_elip']),x['cross_rx_diat']*sin(2*x['cross_rx_rot'])*cos(2*x['cross_rx_elip']),x['cross_rx_diat']*sin(2*x['cross_rx_elip'])))  # 1 x 4
    Dm = np.array((1,x['mol_rx_diat']*cos(2*x['mol_rx_rot'])*cos(2*x['mol_rx_elip']),x['mol_rx_diat']*sin(2*x['mol_rx_rot'])*cos(2*x['mol_rx_elip']),x['mol_rx_diat']*sin(2*x['mol_rx_elip'])))  # 1 x 4
    Dt = np.array((1,x['total_rx_diat']*cos(2*x['total_rx_rot'])*cos(2*x['total_rx_elip']),x['total_rx_diat']*sin(2*x['total_rx_rot'])*cos(2*x['total_rx_elip']),x['total_rx_diat']*sin(2*x['total_rx_elip'])))  # 1 x 4
    Dl = np.array((1,x['lowgain_rx_diat']*cos(2*x['lowgain_rx_rot'])*cos(2*x['lowgain_rx_elip']),x['lowgain_rx_diat']*sin(2*x['lowgain_rx_rot'])*cos(2*x['lowgain_rx_elip']),x['lowgain_rx_diat']*sin(2*x['lowgain_rx_elip'])))
    # Matt: Changed to Matrix Multiplication
    Mir1 = np.dot(RMueller(x['mirror_rot']),np.dot(VWP(0,x['mirror_phase']),np.dot( DiatU(x['mirror_diat']*np.array((1,0,0))), RMueller(-x['mirror_rot']) ) ) )   # 4 x 4
    Mir2 = RevMueller(Mir1)  # 4 x 4
    # Matt: Changed to Matrix Multiplication
    Twin1 = np.dot(RMueller(x['window_diat_rot']),np.dot(DiatU(x['window_diat']*np.array((1,0,0))) ,  \
        np.dot(RMueller(-x['window_diat_rot']),VWP(x['window_phase_rot'],x['window_phase']) ) ) ) # 4 x 4
    Twin2 = RevMueller(Twin1)  # 4 x 4

    for ai in range(len(thQWP)):
        QWP1 = VWP(-thQWP[ai]+x['qwp_rot_0'],
                   x['qwp_phase_0']+x['qwp_phase_2']*np.cos(-2*(thQWP[ai]-x['qwp_rot_0']))+x['qwp_phase_1']*
                   np.cos(-thQWP[ai]+x['qwp_rot_1']))  # 4 x 4
        QWP2 = RevMueller(QWP1)  # 4 x 4
        Si[:,ai] = np.dot(np.dot(np.dot(Twin1,Mir1),QWP1),Stx)
        DxV[ai,:] = np.dot(np.dot(np.dot(Dx,QWP2),Mir2),Twin2)
        DtV[ai,:] = np.dot(np.dot(np.dot(Dt,QWP2),Mir2),Twin2)
        DmV[ai,:] = np.dot(np.dot(np.dot(Dm,QWP2),Mir2),Twin2)

    # convert to T_Array, so we can process properly
    Si = hau.T_Array(Si)
    DxV = hau.T_Array(DxV)
    DtV = hau.T_Array(DtV)
    DmV = hau.T_Array(DmV)
    out = { 'Si': Si,  'Dm': Dm, 'Dx': DxV, 'Dt': DtV, 'DmV': DmV, \
            'QWP1' : QWP1, 'QWP2': QWP2, 'Mir1': Mir1, 'Mir2': Mir2,
            'thQWP': thQWP, 'Twin1' : Twin1, 'Twin2' : Twin2, 'Dt': Dt,
            'Stx': Stx }

#    return out
    return [Si,DxV,DtV,DmV]

def RMueller(theta):
    """ theta - 4x4 array 
        Outputs the 4x4 Mueller matrix for a linear rotation of angle theta (in radians)
    """

    outMat =  np.array( ((1,0, 0, 0), 
                         (0, cos(2*theta), -sin(2*theta), 0),
                         (0,sin(2*theta),cos(2*theta),0), 
                         (0, 0, 0, 1)))  # 4 x 4
    return outMat

def VWP(theta,Gamma):
    """ 
    Outputs the 4x4 Mueller matrix of a linear wave plate of orientation theta (in radians) and 
    phase shift Gamma (in radians).

    2012/01/13 Updated 34 and 43 sign on phase shift matrix

    note: matrix multiply of 3  4x4 arrays
    """
    outMat = np.dot(np.dot (np.array(( (1, 0, 0, 0), 
                                       ( 0, np.cos(2*theta), -np.sin(2*theta), 0), 
                                       (0, np.sin(2*theta), np.cos(2*theta), 0), 
                                       (0, 0, 0,1))),  
                            np.array(( (1,0, 0,0),
                                       (0, 1, 0, 0), 
                                       (0, 0, np.cos(Gamma), np.sin(Gamma)), 
                                       (0, 0, -np.sin(Gamma), np.cos(Gamma)) ))),
                    np.array(( (1, 0, 0, 0),
                               (0, np.cos(-2*theta), -np.sin(-2*theta), 0), 
                               (0, np.sin(-2*theta), np.cos(-2*theta), 0),
                               ( 0, 0, 0, 1))))  # 4 x 4
    return outMat

def DiatU(D):
    """ % Outputs the 4x4 unnormalized Mueller matrix of a diattenuator having a 3 x1 diattenuation vector D. 

    (TODO: transpose may not be needed for 3x1 vector )
    """
    # 
    Dm = np.linalg.norm(D)
    if Dm != 0:
        MD = np.zeros( (4,4) )
        MD[0,1:] = np.transpose(D)
        MD[:,0] = np.append(1, D)
        MD += np.dot(np.sqrt(1-Dm**2),
                     np.array( ( (0 ,0, 0, 0),
                                 (0, 1, 0, 0), 
                                 (0, 0, 1, 0), 
                                 (0, 0, 0, 1)))) + \
            np.dot(
                np.transpose(np.array([np.dot( (1-np.sqrt(1-Dm**2)),np.append(0,D))])), 
                np.array([np.dot(np.append(0, D),1/Dm**2)]) )   
            # Matt: added np.array after first dot product.  
            # Output was otherwise a 1D array and result was a scalar instead of 4x4
    else:           
        MD = np.identity(4)  # 4 x 4
    return MD

def RevMueller(Mi):
    """
    Mr = RevMueller(Mi)
# For a system described by the Mueller matrix Mi propagating in the k
# direction, this function returns the Mueller matrix, Mr for propagating
# through the same system in the -k direction.
"""
    Mr = np.transpose(Mi)* \
                np.array(( (1, 1, -1, 1), (1, 1, -1, 1), (-1, -1, 1, -1),(1, 1, -1, 1)))  
    # Matt: Changed to element wise multiplication
    return Mr



def HSRL_FormMeasMatrix(Si,Dr,varCount):
    """
    Computes the Measurement Matrix measM for a given set of transmitted
    Stokes Vector Si (4xN) and received Diattenuation Vectors Dr (Nx4).
    vars indicates the number of scattering matrix variables under
    consideration (generally 2 or 10)
    """

    measM = np.zeros( (varCount,) )  # vars 

    szm = Si.shape[1] -1
    if szm > 0:
        # Numerically Integrate the Measurement
        if varCount == 1:
            # Only use for molecular profile at maximum temporal resolution.
            # Assumes depolarization is effectively zero.
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0])+Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 2:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0])+Si[3,:]*np.transpose(Dr[:,3]))/szm
            measM[1] = np.trapz(Si[1,:]*np.transpose(Dr[:,1])-Si[2,:]*np.transpose(Dr[:,2])
                                  -2*Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 6:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0]))/szm
            measM[1] = np.trapz(Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))/szm
            measM[2] = np.trapz(Si[1,:]*np.transpose(Dr[:,1]))/szm
            measM[3] = np.trapz(Si[2,:]*np.transpose(Dr[:,2]))/szm
            measM[4] = np.trapz(Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))/szm
            measM[5] = np.trapz(Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 7:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0]))/szm
            measM[1] = np.trapz(Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))/szm
            measM[2] = np.trapz(Si[3,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,3]))/szm
            measM[3] = np.trapz(Si[1,:]*np.transpose(Dr[:,1]))/szm
            measM[4] = np.trapz(Si[2,:]*np.transpose(Dr[:,2]))/szm
            measM[5] = np.trapz(Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))/szm
            measM[6] = np.trapz(Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 9:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0]))/szm
            measM[1] = np.trapz(Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))/szm
            measM[2] = np.trapz(Si[2,:]*np.transpose(Dr[:,0])-Si[0,:]*np.transpose(Dr[:,2]))/szm
            measM[3] = np.trapz(Si[1,:]*np.transpose(Dr[:,1]))/szm
            measM[4] = np.trapz(Si[2,:]*np.transpose(Dr[:,1])-Si[1,:]*np.transpose(Dr[:,2]))/szm
            measM[5] = np.trapz(Si[3,:]*np.transpose(Dr[:,1])+Si[1,:]*np.transpose(Dr[:,3]))/szm
            measM[6] = np.trapz(Si[2,:]*np.transpose(Dr[:,2]))/szm
            measM[7] = np.trapz(Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))/szm
            measM[8] = np.trapz(Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 10:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0]))/szm
            measM[1] = np.trapz(Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))/szm
            measM[2] = np.trapz(Si[2,:]*np.transpose(Dr[:,0])-Si[0,:]*np.transpose(Dr[:,2]))/szm
            measM[9] = np.trapz(Si[3,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,3]))/szm
            measM[3] = np.trapz(Si[1,:]*np.transpose(Dr[:,1]))/szm
            measM[4] = np.trapz(Si[2,:]*np.transpose(Dr[:,1])-Si[1,:]*np.transpose(Dr[:,2]))/szm
            measM[5] = np.trapz(Si[3,:]*np.transpose(Dr[:,1])+Si[1,:]*np.transpose(Dr[:,3]))/szm
            measM[6] = np.trapz(Si[2,:]*np.transpose(Dr[:,2]))/szm
            measM[7] = np.trapz(Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))/szm
            measM[8] = np.trapz(Si[3,:]*np.transpose(Dr[:,3]))/szm
        elif varCount == 16:
            for mi in range(4):
                for ki in range(4):
                    measM[mi*4+(ki+1),1] = trapz(Si[ki,:]*np.transpose(Dr[:,mi]))/szm
    else:
        # Calculate measurement for single case
        if varCount == 1:
            measM[0] = np.trapz(Si[0,:]*np.transpose(Dr[:,0])+Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 2:
            measM[0] = (Si[0,:]*np.transpose(Dr[:,0])+Si[3,:]*np.transpose(Dr[:,3]))
            measM[1] = (Si[1,:]*np.transpose(Dr[:,1])-Si[2,:]*np.transpose(Dr[:,2])-2*Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 6:
            measM[0] = (Si[0,:]*np.transpose(Dr[:,0]))
            measM[1] = (Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))
            measM[2] = (Si[1,:]*np.transpose(Dr[:,1]))
            measM[3] = (Si[2,:]*np.transpose(Dr[:,2]))
            measM[4] = (Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))
            measM[5] = (Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 7:
            measM[0] = (Si[0,:]*np.transpose(Dr[:,0]))
            measM[1] = (Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))
            measM[2] = (Si[3,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,3]))
            measM[3] = (Si[1,:]*np.transpose(Dr[:,1]))
            measM[4] = (Si[2,:]*np.transpose(Dr[:,2]))
            measM[5] = (Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))
            measM[6] = (Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 9:
            measM[0] = (Si[0,:]*np.transpose(Dr[:,0]))
            measM[1] = (Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))
            measM[2] = (Si[2,:]*np.transpose(Dr[:,0])-Si[0,:]*np.transpose(Dr[:,2]))
            measM[3] = (Si[1,:]*np.transpose(Dr[:,1]))
            measM[4] = (Si[2,:]*np.transpose(Dr[:,1])-Si[1,:]*np.transpose(Dr[:,2]))
            measM[5] = (Si[3,:]*np.transpose(Dr[:,1])-Si[1,:]*np.transpose(Dr[:,3]))
            measM[6] = (Si[2,:]*np.transpose(Dr[:,2]))
            measM[7] = (Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))
            measM[8] = (Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 10:
            measM[0] = (Si[0,:]*np.transpose(Dr[:,0]))
            measM[1] = (Si[1,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,1]))
            measM[2] = (Si[2,:]*np.transpose(Dr[:,0])-Si[0,:]*np.transpose(Dr[:,2]))
            measM[9] = (Si[3,:]*np.transpose(Dr[:,0])+Si[0,:]*np.transpose(Dr[:,3]))
            measM[3] = (Si[1,:]*np.transpose(Dr[:,1]))
            measM[4] = (Si[2,:]*np.transpose(Dr[:,1])-Si[1,:]*np.transpose(Dr[:,2]))
            measM[5] = (Si[3,:]*np.transpose(Dr[:,1])+Si[1,:]*np.transpose(Dr[:,3]))
            measM[6] = (Si[2,:]*np.transpose(Dr[:,2]))
            measM[7] = (Si[3,:]*np.transpose(Dr[:,2])-Si[2,:]*np.transpose(Dr[:,3]))
            measM[8] = (Si[3,:]*np.transpose(Dr[:,3]))
        elif varCount == 16:
            for mi in range(4):
                for ki in range(4):
                    measM[(mi)*4+(ki+1),1] = Si[ki,:]*np.transpose(Dr[:,mi])


    return measM
