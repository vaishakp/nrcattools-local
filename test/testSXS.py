''' Test SXS waveform load against LALSimulation '''



import numpy as np
import sys
import os
import h5py

# nr-catalog-tools
cwd = os.getcwd()

libpath = f'{cwd}/../'

if libpath not in sys.path:
    sys.path.append(libpath)

from nrcatalogtools.sxs import SXSCatalog

#import matplotlib.pyplot as plt

# unittest funcs
from helper import *
import unittest

# waveformtools
from waveformtools.waveforms import modes_array
from waveformtools.waveformtools import lengtheq

# pycbc
from pycbc.waveform import td_approximants
from pycbc.types.timeseries import TimeSeries
from pycbc.filter.matchedfilter import match
from pycbc.waveform.utils import coalign_waveforms
from pycbc import pnutils




######################################
# Simulation properties
######################################

# Simulation name
sim_name = 'SXS:BBH:0001'
# Parameters
M = 40
D = 1000
inc = np.pi/6
coa_phase = np.pi/4
delta_t = 1./2048
# Extrinsic parameters:
f_lower = 20
f_lower_at_1MSUN = f_lower/M

# Convention
# hp1, hx1, h1... : nrcat waveforms
# wfa1, hp2, hx2, ... : waveformtools waveforms


#######################################
# Fetch waveform using nr-catalog-tools
#######################################

#sc = sxs.Catalog.load(download=True)
#rc = RITCatalog.load(verbosity=5, download=True)
message('Loading SXS waveform through nrcatalogtools...')
sxs1 = SXSCatalog.load(download=True)
#mc = MayaCatalog.load(verbosity=5)

#mwf = mc.get(sim_name)
sxsw = sxs1.get(sim_name)


hpc = sxsw.get_td_waveform(total_mass=M, distance=D, inclination=inc,
                    coa_phase=coa_phase, delta_t=delta_t
                    )
hpc_pycbc = hpc # mwf.to_pycbc(hpc)
hp1, hx1 = hpc_pycbc.real(), hpc_pycbc.imag()
h1 = hp1 + 1j*hx1

#plt.plot(hp1.sample_times, hp1)
#plt.plot(hx1.sample_times, hx1)
#plt.grid()
#plt.show()



#########################################
# Fetch waveform using waveformtools
#########################################


fdir = "/home/runner/.cache/sxs/SXS:BBH:0001v3/Lev5"
fname = 'rhOverM_Asymptotic_GeometricUnits_CoM.h5'

message('Loading SXS waveform through waveformtools...')

wfa1 = modes_array(label='sxs_001', spin_weight=-2)
wfa1.file_name = fname
wfa1.data_dir = fdir
wfa1.load_modes(ftype='SpEC', var_type='strain', ell_max='auto', resam_type='auto')
wfa1.get_metadata()


taxis, hp2, hx2 = wfa1.to_td_waveform(Mtotal=M, 
                                      distance=D, 
                                      incl_angle=inc, 
                                      delta_t=delta_t)


h2 = hp2 + 1j*hx2



#plt.plot(t, hp2, color=[0,0.7071,1])
#plt.plot(t, hx2, color=[0.1,0,0])
#plt.show()

###################
# Initiate tests
##################

# Equalize lengths

h = lengtheq(np.array(h1), h2)
_, h2, flag = h
print(f'Waveform {flag} length changed')
hp2 = h2.real
hx2 = h2.imag

hp2_ts = TimeSeries(hp2, delta_t=delta_t)
hx2_ts = TimeSeries(hx2, delta_t=delta_t)

mp, sp = match(hp1, hp2_ts)
mx, sx = match(hx1, hx2_ts)

wf1_p, wf2_p = coalign_waveforms(hp1, hp2_ts)
wf1_x, wf2_x = coalign_waveforms(hx1, hx2_ts)

# Normalize the arrays
wf1 = np.array(wf1_p) + 1j*np.array(wf1_x)
wf2 = np.array(wf2_p) + 1j*np.array(wf2_x)

n1 = np.sqrt(np.dot(wf1, np.conjugate(wf1)))
n2 = np.sqrt(np.dot(wf2, np.conjugate(wf2)))

wf1 = wf1/n1
wf2 = wf2/n2

wf1_p = wf1.real
wf1_x = wf1.imag

wf2_p = wf2.real
wf2_x = wf2.imag

wf1_p = TimeSeries(wf1_p, delta_t)
wf1_x = TimeSeries(wf1_x, delta_t)
wf2_p = TimeSeries(wf2_p, delta_t)
wf2_x = TimeSeries(wf2_x, delta_t)



class TestSXS(unittest.TestCase):
    ''' Test loading of SXS waveforms '''
    
    def test_waveforms(self):
        ''' Test the SXS loading of waveforms against 
        that loading using lalsimulation. Tested are RMS errors, maximum deviation and mismatches'''
        
        # L2 errors
        Res_p, Amin_p, Amax_p = RMSerrs(np.array(wf1_p), np.array(wf2_p))
        Res_x, Amin_x, Amax_x = RMSerrs(np.array(wf1_x), np.array(wf2_x))
        
        #Amin_p/=A1max
        #Amin
        # Match
        match_p, shift_p = match(wf1_p, wf2_p)
        match_x, shift_x = match(wf1_x, wf2_x)

        mismatch_p = 100*(1-match_p)
        mismatch_x = 100*(1-match_x)
        
        prec = 1
        # RMS error should be less than 0.1 x Amax(wf1)
        self.assertAlmostEqual(Res_p, 0, prec, f"The RMS error between the + components of the waveforms must be atmost 0.1 times Max amplitude of the normalized waveform")
        self.assertAlmostEqual(Res_x, 0, prec, f"The RMS error between the x components of the waveforms must be almost 0.1 times Max amplitude of the normalized waveform")
        
        prec = 0
        # Max relative point-wise deviation w.r.t Amax(wf1) should be less than 1 (100)%
        self.assertAlmostEqual(np.absolute(Amin_p), 0, prec, f"The maximum lower deviation between the + components of the waveforms must be almost 100%")
        self.assertAlmostEqual(np.absolute(Amax_p), 0, prec, f"The maximum upper deviation between the x components of the waveforms must be almost 100%")
        
        self.assertAlmostEqual(np.absolute(Amax_p), 0, prec, f"The maximum upper deviation between the + components of the waveforms must be almost 0")
        self.assertAlmostEqual(np.absolute(Amax_p), 0, prec, f"The maximum upper deviation between the x components of the waveforms must be almost 0")
        
        prec = 1
        # Mismatch should be less than 0.1%
        self.assertAlmostEqual(mismatch_p, 0, prec, f"The mismatch between the + components of the waveforms must be almost 0.1%")
        self.assertAlmostEqual(mismatch_p, 0, prec, f"The mismatch between the x components of the waveforms must be almost 0.1%")

       
        prec=1
        # Full array
        np.testing.assert_almost_equal(wf1, wf2, prec)

if __name__ == '__main__':
    unittest.main()
