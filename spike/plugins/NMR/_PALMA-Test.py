#!/usr/bin/env python
# coding: utf-8

import platform
import ctypes
import os.path as op
from datetime import datetime
import time
import multiprocessing as mp
import warnings

import spike
from spike.File.BrukerNMR import Import_1D
from spike.plugins.NMR.PALMA import Import_DOSY
from spike.Tests import filename, directory

print('Run date:', datetime.now().isoformat() )

# start chronometer
t0 = time.time()
# ### Import dataset
selected = filename("HEKmedia/30/ser")

d2 = Import_DOSY(selected)
d2.axis1.dmin = 10        # diffusion axis limits, values are in µm²/cm
d2.axis1.dmax = 10000     # only indicative à this level
d2.filename = selected
d2.pulprog = d2.params['acqu']['$PULPROG']
print(f"Gradients values from {d2.axis1.qvalues[0]} to {d2.axis1.qvalues[-1]} G/cm in {d2.size1} steps")
print (d2.params['acqu']['title'])

D2 = d2.copy() # copy the imported data-set to another object for processing
# bk_ftF2 and bk_ftF1 (define in the Bruker plugin) find which FT to apply depending on FnMODE
D2.apod_sin(maxi=0.1,axis='F2').zf(1,2).bk_ftF2().bk_pk()  # chaining  apodisation - zerofill - FT
D2.phase(-64.3,  37.0).real()     # phase it, then throw im part
D2.axis2.currentunit = 'ppm'
#D2.display(scale="auto", autoscalethresh=6.0, title="%s %s"%(selected,d2.pulprog))  # chain  set to ppm unit - and display

# simple linear baseline correction using extreme points
for i in range(D2.size1):
    r=D2.row(i)
    r.bcorr(method='spline', xpoints=[0.106, 9.838], nsmooth=5, xpunits='current',)
    D2.set_row(i,r)



# ## PALMA processing

############################################## First set your parameters
# 
# Processing time is proportional to $N$ `x nbiter x finalsize`
# 
# where $N$ is the number of processed columns (chosen with `miniSNR`)

# Diffusion axis
finalsize = 256  # The final of the DOSY on the diffusion axis
Dmin = 30        # the minimum diffusion coefficient (in µm2/sec) typical is 1-100
Dmax = 10000     # the maximum diffusion coefficient (in µm2/sec) typical is 5000-50000

# Processing
nbiter=1000       # number of iterations - 5000-20000 - the more the better (usually)
lamda=0.01        # weight between "pure" MaxEnt (1.0) and "pure" l1 (0.0), 0.01 to 0.1 are "good values"

# Optionnal parameters
miniSNR = 32          # minimum SNR in column to do processing - 32 is optimal - do not go below 8

# uncertainty=1.2  # if >1 allows more room for algo when data quality is poor
# precision=1e-08  # stopping criterium

# MultiProcessing
# the processing can be lengthy,  so use can use parralelize the program
# if you do not want to use the mp capability, set NProc to 1
NProc = 4         # here for 4 cores - adapt to your own requirements
#############################################################################

if NProc > 1:
    # this forces MKL -if there- to use only one thread (usually faster ! and probably required if you intend using multiprocessing)
    dllnames = dict(Linux='libmkl_rt.so', Darwin='libmkl_rt.dylib', Windows='libmkl_rt.dll')
    dllname = dllnames[platform.system()]

    mkl_rt = ctypes.CDLL(dllname)
    mkl_rt.MKL_Set_Num_Threads(1)
    print ("MKL stopped")

    mppool = mp.Pool(processes=NProc)

else:
    mppool = None

D2.prepare_palma(finalsize, Dmin, Dmax)      # set-up diffusion axis

with warnings.catch_warnings():
    warnings.simplefilter(action="ignore", category="RuntimeWarning", lineno=87)
# action:message:category:module:line
    DD2 = D2.do_palma(nbiter=nbiter, miniSNR=miniSNR, lamda=lamda, mppool=mppool)  # launch computation

if mppool: mppool.close()  # clean-up

# ### And display result
# Diffusion coefficients are along the vertical axis. Values are in $\mu m^2 / sec$ which is $10^{-12} m^2 / sec$.

DD2.axis1.currentunit = 'Diff'

#DD2.display(scale="auto", autoscalethresh=6.0, title=f"{selected} PALMA DOSY processing))



# ## Save the data-set
#DD2.save('processed_Dosy.gs2')

proc_time =  round(time.time()-t0)
print (f"The whole processing took {proc_time} sec : {round(proc_time/60,2)} minutes")

import unittest
class FitTests(unittest.TestCase):
    def test_DOSY(self):
        self.assertAlmostEqual(D2.axis1.dfactor, 1162859.287848, 3)
        #Test results
        c1 = DD2.colc(0.8420)  # extract a single component column at 0.8429 ppm
        print(c1[int(c1.axis1.dtoi(2685))], c1[int(c1.axis1.dtoi(4500))], c1[int(c1.axis1.dtoi(1500))])
        # c1[int(c1.axis1.dtoi(2685))], c1[int(c1.axis1.dtoi(4500))], c1[int(c1.axis1.dtoi(1500))]
        # (424860.47585149534, -213.78871091846162, -280.7466564428571)   
        self.assertTrue(c1[int(c1.axis1.dtoi(2685))] > 400000)   # coordinates in mu^2/sec
        self.assertTrue(c1[int(c1.axis1.dtoi(4500))] < 4000)
        self.assertTrue(c1[int(c1.axis1.dtoi(1500))] < 4000)

        c2 = DD2.colc(3.5958)  # then a polydisperse component
        print(c2[int(c2.axis1.dtoi(3017))],c2[int(c2.axis1.dtoi(1400))], c2[int(c2.axis1.dtoi(384))])
        self.assertTrue(c2[int(c2.axis1.dtoi(3017))] > 100000)
        self.assertTrue(c2[int(c2.axis1.dtoi(1386))] < 40000)
        self.assertTrue(c2[int(c2.axis1.dtoi(384))] > 120000)
        # self.assertTrue(c2[int(c2.axis1.dtoi(384))] < 200000)
        # c2[int(c2.axis1.dtoi(3017))],c2[int(c2.axis1.dtoi(1386))], c2[int(c2.axis1.dtoi(384))]
        # (147903.9573949227, 22981.56238030312, 156769.0940550826)
        # 536015.7981013925 -182.74007087275956 -241.24269380861224
        # 195166.05636473492 28147.27290978278 211751.14251187275


if __name__ == '__main__':
    unittest.main()
