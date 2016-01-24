#!/usr/bin/env python 
# encoding: utf-8

"""
This a test suite for NPK.
  
Note that all functionalities  are not tested yet !
"""
from __future__ import print_function
import unittest
import tempfile
import math
from v1.Kore import *

error=0
#compatibility(locals())
class TextNPKv1(unittest.TestCase):
    def setUp(self):
        self.verbose = 1
    def announce(self):
        if self.verbose >0:
            print("\n========",self.shortDescription(),'===============')
    def test_benchmark(self):
        """
     benchmark for gifa
     verbose 1    # check what you're doing
     benchmark for NPK
       run on a 512 x 2k real data set.
       Processing performed is :
       - Exponen. Broading in F2
       - Real FT, phasing, in F22
       - 5 points Spline base-line correction in F2
       - cosine apodization in F1, zero-filling in F1
       - Real FT and phasing in F1.
     The displayed spectrum at the end is 512 by 1k real points
     results are (some figures date from version 1.0 !) :
        """

    # first prepare a pseudo-data for benchmarking GIFA
        import time
        dim(2)
        chsize(128, 128)
        one()
        sin(0.0, "f12")
        revf("f1")
        revf("f2")  # put something in it
        chsize(512, 2048)
        itype(0)                         # and zero-fill
        put("data")

        t=time.clock()
        for i in range(10):
          get("data")

          em(0.0, 1.0)
          revf("f2")
          rft("f2")                                    # main "f2" processing
          phase(20.0, 20.0, "f2")
          real("f2")               #throwing Imaginary parts is needed here
          bcorr(2, 2, "f2", [10, 100, 300, 600, 900])                       # for baseline corr.
          sin(0.0, "f1")
          chsize(get_si1_2d()*2,get_si2_2d())
          rft("f1")                              # main "f1" proc
          phase(20.0, 20.0, "f1")       #done
        print('Total TIME:',time.clock()-t)

        print(' ')
        print('Typical results are :')
        print('=====================')
        print(' 7.00"   - G4 1.68 MHz - Mac Os 10.4 (PowerBook)')
        print(' 3.00"   - macintel mono core - Mac Os 10.4 (MacMini)')
        print(' 2.80"   - PIV 2.8 GHz hyperthreading - Linux')
        print(' 2.45"   - AMD 64  Opteron 2.2 GHz - Linux')
        print(' 1.52"   - Core 2 duo 2.8 GHz - Mac Os 10.5 (powerbook)')
        print(' 0.81"   - Xeon 5570 2.66GHz - Mac Os 10.6 (MacPro)')

    def test_basic(self):
        """test the basic NPK mechanisms
        """
        dim(2)
        ok = get_dim() == 2
        dim(1)
        ok = ok and get_dim() == 1
        self.report_result(ok, "Basic context settings")

    def test_cache(self):
        """test the cache sub-system
        """
        self.report_result(1,"dedans")
        t11 = 1000
        t12 = 248  
        t22 = 1000
        t13 = 40   
        t23 = 120   
        t33 = 80

        # first test if there is enough room
        # total = ( 4*2*power2(t11-1)+4*1024)                                         #1d
        # total = ( total + 4*4*power2(t12-1)*power2(t22-1) +4096 )                   #2d
        # total = ( total + 4*8*power2(t13-1)*power2(t23-1)*power2(t33-1) +4096)      #3d
        # total = ( total + 4*power2(t13-1)*power2(t23-1)*power2(t33-1) +4096)        #newfilec
        # total = (total /1024)   # in ko
        total = ( 4*2*pow(t11-1,2)+4*1024)                                         #1d
        total = ( total + 4*4*pow(t12-1,2)*pow(t22-1,2) +4096 )                   #2d
        total = ( total + 4*8*pow(t13-1,2)*pow(t23-1,2)*pow(t33-1,2) +4096)      #3d
        total = ( total + 4*pow(t13-1,2)*pow(t23-1,2)*pow(t33-1,2) +4096)        #newfilec
        total = (total /1024)   # in ko
        # this will be precise only if txx values are slightly smaller than power of 2

    #    print "Be carefull that this test will create large temporary files on ",tempfile.gettempdir()
    #    print "for a total of about",int(total), "ko"
    #    print "Please make sure that this space is available !"

        # create dummy file_names
        file1 = tempfile.mktemp('.gifatemp')
        file2 = tempfile.mktemp('.gifatemp')
        file3 = tempfile.mktemp('.gifatemp')
        file4 = tempfile.mktemp('.gifatemp')
        file5 = tempfile.mktemp('.gifatemp')
        file6 = tempfile.mktemp('.gifatemp')
        file7 = tempfile.mktemp('.gifatemp')
        file8 = tempfile.mktemp('.gifatemp')

        #
        dim(1)
        chsize(t11)
        one()
        mult(100)
        sin(0)
        writec(file1)
        dim(2)
        chsize(t12, t22)
        one()
        mult(1000)
        sin(0, "f12")
        writec(file2)
        dim(3)
        chsize(t13,t23,t33)
        one()
        mult( 10000)
        sin( 0, "f123")
        writec(file3)

        self.report_result(1,"Write")

        #======read test====
        #
        # read back
        dim(1); chsize(512); zero()
        read(file1)
        com_max()
        ok = (get_si1_1d() == t11) and (math.fabs(geta_max(1)-100)<1) and (math.fabs(geta_max(2)) < 1)
        print(ok)
        dim(2); chsize( 100, 100); zero()
        read(file2)
        com_max()
        ok = (ok and (get_si1_2d()*get_si2_2d() == (t12*t22)) and (math.fabs(geta_max(1))-1000)<1 and math.fabs(geta_max(2)) < 1)
        print(ok)
        dim(3); chsize( 10, 10, 10); zero()
        read(file3)
        com_max()
        ok = (ok and (get_si1_3d()*get_si2_3d()*get_si3_3d()) == (t13*t23*t33))
        print(ok)
        ok = (ok and math.fabs(geta_max(1)-10000)<1 and math.fabs(geta_max(2)) < 1)
        print(geta_max(1))
        print(geta_max(2))
        self.report_result(ok, "Read")

        join(file3)
        join(file1)
        join(file2)
        print("listfilec() yet to be done")
    #    dataset()

       # report_result( (get_c_sizef1()*get_c_sizef2() == t12*t22), "Join")

        pi_2 = (2*math.atan(1))
        # sqrt(2)/2 = cos(pi_2 / 2) is equal to 0.7071 ! , its square is 0.5 !!
        # result are tested at 1.4%, this is because the number of points

        # ==========test getc 1d-1d
        join(file1)
        td = (2*int(t11/4) +1) ; tf = (t11-10)
        test = (100*math.cos(pi_2*tf/t11))
        dim(1); getc( td, tf)
        ok = (get_si1_1d() == ( tf - td +1))
        com_max()
        ok = (ok and math.fabs(geta_max(1)-70.71) < 1 and math.fabs(geta_max(2)-test)<1 )
        self.report_result(ok, "Getc 1d - 1d")
        okok = ok

        # ==========test getc 2d-1d
        join( file2)
        td = 1;     tf = t12
        test = (1000*math.cos(pi_2*tf/t12))
        dim( 1); getc( "f1", 1,  td, tf)
        ok = (get_si1_1d() == ( tf - td +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-1000) < 1 and math.fabs(geta_max(2)-test)< 1 )

        td = (2*int(t22/4) +1);     tf = (t22-10)
        test = (707.1*math.cos(pi_2*tf/t22))
        getc( "f2", (t12/2), td, tf)
        ok = (ok and get_si1_1d() == ( tf - td +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-500) < 10 and math.fabs(geta_max(2)-test)<10 )

        report_result( ok, "Getc 2d - 1d")
        okok = (okok and ok)

        # ==========test getc 3d-1d
        join( file3)
        td = 1;     tf = t13
        test = (10000*math.cos(pi_2*tf/t13))
        dim( 1); getc ("f1", 1, 1, td, tf)
        ok = (get_si1_1d() == ( tf - td +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-10000) < 1 and math.fabs(geta_max(2)-test)< 1 )

        td = (2*int(t23/4) +1);     tf = (t23-10)
        test = (7071*math.cos(pi_2*tf/t23))
        getc( "f2", (t13/2), 1, td, tf)
        ok = (ok and get_si1_1d() == ( tf - td +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-5000) < 100 and math.fabs(geta_max(2)-test)<100 )

        td = (2*int(t33/4) +1);     tf = (t33-10)
        test = (5000*math.cos(pi_2*tf/t33))
        getc( "f3", (t13/2), (t23/2), td, tf)
        ok = (ok and get_si1_1d() == ( tf - td +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-5000*0.7071) < 100 and math.fabs(geta_max(2)-test)<100 )

        self.report_result( ok, "Getc 3d - 1d")
        okok = (okok and ok)

        # ==========test getc 2d-2d
        join( file2)
        td1 = (2*int(t12/4) +1);     tf1 = (t12-10)
        td2 = (2*int(t22/4) +1);     tf2 = (t22-10)
        test = (1000*math.cos(pi_2*tf1/t12)*math.cos(pi_2*tf2/t22))
        dim( 2); getc( td1, td2, tf1, tf2)
        ok = (get_si1_2d() == ( tf1 - td1 +1) and get_si2_2d() == ( tf2 - td2 +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-500) < 10 and math.fabs(geta_max(2)-test)< 10 )

        self.report_result( ok, "Getc 2d - 2d")
        okok = (okok and ok)

        # ==========test getc 3d-2d
        join( file3)
        td1 = 1;     tf1 = t23
        td2 = 1;     tf2 = t33
        test = (10000*math.cos(pi_2*tf1/t23)*math.cos(pi_2*tf2/t33))
        dim( 2); getc( "f1", 1, td1, td2, tf1, tf2)
        ok = (get_si1_2d() == ( tf1 - td1 +1) and get_si2_2d() == ( tf2 - td2 +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-10000) < 1 and math.fabs(geta_max(2)-test)< 1 )

        td1 = (2*int(t13/4) +1);     tf1 = (t13-10)
        td2 = (2*int(t33/4) +1);     tf2 = (t33-10)
        test = (10000*math.cos(pi_2*tf1/t13)*math.cos(pi_2*tf2/t33))
        dim( 2); getc( "f2", 1, td1, td2, tf1, tf2)
        ok = (ok and get_si1_2d() == ( tf1 - td1 +1) and get_si2_2d() == ( tf2 - td2 +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-5000) < 200 and math.fabs(geta_max(2)-test)< 100 )

        td1 = (2*int(t13/4) +1);     tf1 = (t13-10)
        td2 = (2*int(t23/4) +1);     tf2 = (t23-10)
        test = (7071*math.cos(pi_2*tf1/t13)*math.cos(pi_2*tf2/t23))
        dim( 2); getc( "f3", (t33/2), td1, td2, tf1, tf2)
        ok = (ok and get_si1_2d() == ( tf1 - td1 +1) and get_si2_2d() == ( tf2 - td2 +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-5000*0.7071) < 100 and math.fabs(geta_max(2)-test)< 100 )

        self.report_result( ok, "Getc 3d - 2d")
        okok = (okok and ok)

        # ==========test getc 3d-3d
        join( file3)
        td1 = (2*int(t13/4) +1);     tf1 = (t13-10)
        td2 = (2*int(t23/4) +1);     tf2 = (t23-10)
        td3 = (2*int(t33/4) +1);     tf3 = (t33-10)
        test = (10000*math.cos(pi_2*tf1/t13)*math.cos(pi_2*tf2/t23)*math.cos(pi_2*tf3/t33))
        dim (3); getc( td1, td2, td3, tf1, tf2, tf3)
        ok = (get_si1_3d() == ( tf1 - td1 +1) and get_si2_3d() == ( tf2 - td2 +1) and get_si3_3d() == ( tf3 - td3 +1))
        com_max() 
        ok = (ok and math.fabs(geta_max(1)-5000*0.7071) < 200 and math.fabs(geta_max(2)-test)< 100 )
        report_result( ok, "Getc 3d - 3d")
        okok = (okok and ok)

        self.report_result( okok, "Getc")


        # make zero files to put in
        dim( 1); chsize( t11);		zero(); writec( file4)
        dim( 2); chsize( t12, t22);		zero(); writec( file5)
        dim( 3); chsize( t13, t23, t33);	zero(); writec( file6)

        #=========================
        # ==========test putc 1d
        join( file4)
        td = 1;     tf = (t11/2)
        dim(1); chsize( tf); one()
        putc( td, tf)
        zero()
        getc(td,t11)
        com_max()
        ok = (geta_max(1) == 1.0 and math.fabs(geta_max(2)) < 1.0e-6)
        getc( (tf+1), t11)
        com_max()
        ok = (ok and math.fabs(geta_max(1)) < 1.0e-6 and math.fabs(geta_max(2)) < 1.0e-6)
        self.report_result( ok, "Putc 1d - 1d")
        okok = ok

        # ==========test putc 2d-1d
        join( file5)
        td = 1;     tf = (t12/2)
        dim(1); chsize(tf); one()
        putc( "f1", 1 , td, tf)
        zero()
        getc( "f1", 1, td, t12)
        com_max()
        ok = (geta_max(1) == 1.0 and math.fabs(geta_max(2) < 1.0e-6))
        getc( "f1", 1, (tf+1), t12)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        td = 1;     tf = (t22/2)
        dim(1); chsize( tf); one()
        putc( "f2", (t12/2 +1),  td, tf)
        zero()
        getc( "f2", (t12/2 +1), td, t22)
        com_max()
        ok = (ok and geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f2", (t12/2 +1), (tf+1), t22)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        self.report_result( ok, "Putc 2d - 1d")
        okok = (okok and ok)


        # ==========test putc 3d-1d
        join( file6)
        td = 1;     tf = (t13/2)
        dim(1); chsize(tf); one()
        putc( "f1", 1, 1, td, tf)
        zero()
        getc( "f1", 1, 1,td,t13)
        com_max()
        ok = (geta_max(1) == 1.0 and math.fabs(geta_max(2)) < 1.0e-6)
        getc("f1", 1, 1, (tf+1),t13)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        td = 1;     tf = (t23/2)
        dim(1); chsize( tf); one()
        putc( "f2", 1, (t13/2+1), td, tf)
        zero()
        getc( "f2", 1, (t13/2+1),td, t23)
        com_max()
        ok = (ok and geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f2", 1, (t13/2+1), (tf+1), t23)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        td = 1;     tf = (t33/2)
        dim(1); chsize( tf); one()
        putc( "f3", (t13/2+1), (t23/2+1), td, tf)
        zero()
        getc( "f3", (t13/2+1), (t23/2+1),td,t33)
        com_max()
        ok = (ok and geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f3", (t13/2+1),  (t23/2+1), (tf+1),t33)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        self.report_result( ok, "Putc 3d - 1d")
        okok = (okok and ok)

        # ==========test putc 2d-2d
        join( file5)
        td1 = 1;     tf1 = (t12/2)
        td2 = 1;     tf2 = (t22/2)
        dim(2); chsize(tf1, tf2); one()
        putc( td1, td2, tf1, tf2)
        zero()
        getc( td1, td2, t12, t22)
        com_max()
        ok = (geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( (tf1+1), (tf2+1), t12,t22)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        self.report_result( ok, "Putc 2d - 2d")
        okok = (okok and ok)

        # ==========test putc 3d-2d
        join( file6)
        td1 = 1;     tf1 = (t23/2)
        td2 = 1;     tf2 = (t33/2)
        dim(2); chsize(tf1, tf2); one()
        putc( "f1", 1, td1, td2, tf1, tf2)
        zero()
        getc( "f1", 1, td1, td2, t23, t33)
        com_max()
        ok = (geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f1", 1, (tf1+1) ,(tf2+1),t23,t33)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        td1 = 1;     tf1 = (t13/2)
        td2 = 1;     tf2 = (t33/2)
        dim(2); chsize(tf1, tf2); one()
        putc( "f2", (t23/2+1), td1, td2, tf1, tf2)
        zero()
        getc( "f2", (t23/2+1), td1, td2, t13, t33)
        com_max()
        ok = (ok and geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f2", (t23/2+1), (tf1+1), (tf2+1),t13,t33)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        td1 = 1;     tf1 = (t13/2)
        td2 = 1;     tf2 = (t23/2)
        dim(2); chsize(tf1, tf2); one()
        putc( "f3", (t33/2+1), td1, td2, tf1, tf2)
        zero()
        getc( "f3", (t33/2+1), td1, td2, t13, t23)
        com_max()
        ok = (ok and geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( "f3", (t33/2+1), (tf1+1), (tf2+1), t13,t23)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        self.report_result( ok, "Putc 3d - 2d")
        okok = (okok and ok)

        # ==========test putc 3d-3d
        join( file6)
        td1 = 1;     tf1 = (t13/2)
        td2 = 1;     tf2 = (t23/2)
        td3 = 1;     tf3 = (t33/2)
        dim(3); chsize(tf1, tf2, tf3); one()
        putc( td1, td2, td3, tf1, tf2, tf3)
        zero()
        getc( td1, td2, td3, t13, t23, t33)
        com_max()
        ok = (geta_max(1) == 1.0 and geta_max(2) < 1.0e-6)
        getc( (tf1+1), (tf2+1), (tf3+1), t13,t23,t33)
        com_max()
        ok = (ok and geta_max(1) < 1.0e-6 and geta_max(2) < 1.0e-6)

        self.report_result( ok, "Putc 3d - 3d")
        okok = (okok and ok)

        report_result( okok, "Putc")

        #=========================
        # ======== newfilec et proc2d
    
    #    ( file feq itype (tailles) (offset) (sw) (freq) )
        newfilec( file7, 10.1, 0, 
            (t13/4, t23/4, t33/4),
            (0.1,0.2,0.3),
            (1.1,1.2,1.3),
            (2.1,2.2,2.3))
        disjoin()
        join( file7)
    #    filec_status()
        ok = (get_c_dim() == 3 and math.fabs(get_c_offsf1()-0.1)<1.0e-6 and math.fabs(get_c_offsf2()-0.2)<1.0e-6 and math.fabs(get_c_offsf3()-0.3)<1.0e-6)
        ok = (ok and math.fabs(get_c_specwf1()-1.1)<1.0e-6 and math.fabs(get_c_specwf2()-1.2)<1.0e-6 and math.fabs(get_c_specwf3()-1.3)<1.0e-6 )
        ok = (ok and math.fabs(get_c_freq1()-2.1)<1.0e-6 and math.fabs(get_c_freq2()-2.2)<1.0e-6 and math.fabs(get_c_freq3()-2.3)<1.0e-6)
        ok = (ok and get_c_sizef1()==(t13/4) and get_c_sizef2()==(t23/4) and get_c_sizef3()==(t33/4))
        ok = (ok and get_c_type() == 0 and math.fabs(get_c_freq()-10.1)<1.0e-6 )
        dim(1); chsize(t23/4)
        for i in range(1,1+t13/4):
            for j in range(1,1+t33/4):
                one(); mult( i+j); putc( "f2", i, j, 1, (t23/4))
        disjoin()
        join( file7)
    #    filec_status()
        dim(2)
        getc( "f2", (t23/8), 1, 1, (t13/4), (t33/4))
        com_max()
        ok = (geta_max(1) == (t13/4 + t33/4) and geta_max(2) == 2 )

        disjoin()

        self.report_result( ok, "newfilec")
    
        #======= disjoin
        join( file4); disjoin()
        join( file5); disjoin()
        join( file6); disjoin()
        join( file1); disjoin()
        join( file2); disjoin()
        join( file3); disjoin()

        self.report_result (not get_c_joined(), "Disjoin")


        #============= proc2d ====
        print("proc3d yet to be done")
        proc3d(file7, file8, "f1", "chsize(get_si1_2d()/2,get_si2_2d()); mult(-1)",globals() )
        join(file8)
        # procesing results in reducing size in F2 by two
        #filec_status()
        ok = (get_c_sizef1()==(t13/4) and get_c_sizef2()==(t23/8) and get_c_sizef3()==(t33/4))
        getc( "f2", (t23/8), 1, 1, (t13/4), (t33/4))
        com_max()
        ok = (ok and geta_max(2) == -(t13/4 + t33/4) and geta_max(1) == -2 )
        disjoin()

        self.report_result( ok, 'proc3d')

        # ================ Getheader, Putheader

        join( file6) # was file7
        print(file6)
        putheader( "Type", repr(12))
        ok = (getheader("Dim") == "3")
        getheader("Type")
        ok = (ok and getheader("Type") == "12")
        putheader("Coucou", "You can even put string in there")
        print(getheader("Coucou"))
        ok = (ok and getheader("Coucou") ==  "You can even put string in there")
        putheader( "2nd example", "Parameters can be anything")
        ok = (ok and getheader( "2nd example") ==  "Parameters can be anything")
        disjoin()

        join( file6)   # was file7
        ok = (ok and getheader( "Type") == "12")
        disjoin()

        self.report_result( ok, "Getheader, Putheader")

        # ================ addc
        print("addc() yet to be done")
        """
        dim(1); chsize(t11);		one(); writec( file4)
        one(); mult(2)
        addc(file4)
        com_max()
        ok = (geta_max(1) == 3.0)

        dim(2); chsize(t12, t22);		one(); mult( 10);  writec( file5)
        one(); mult( 4)
        addc( file5)
        com_max()
        ok = (ok and geta_max(1) == 14.0)

        dim( 3); chsize(t13, t23, t33);	one(); mult( 100); writec( file6)
        one() ; mult( 6)
        addc( file6)
        com_max()
        ok = (ok and geta_max(1) == 106.0)

        report_result( ok, Addc)
        """

    def test_dataset(self):
        if (get_c_joined() == 0) :
          print("No currently JOINed file")
        else:
            print("Currently JOINed file :")
            print("Filename       = " +repr(get_c_name()))
            print("Dim            = " +repr(get_c_dim()))
            print("Size f1        = " +repr(get_c_sizef1()))
            print("Size f2        = " +repr(get_c_sizef2()))
            print("Size F3        = " +repr(get_c_sizef3()))
            print("1H Frequency   = " +repr(get_c_freq()))
            print("Freq f1        = " +repr(get_c_freq1()))
            print("Freq f2        = " +repr(get_c_freq2()))
            print("Freq F3        = " +repr(get_c_freq3()))
            print("Specw f1       = " +repr(get_c_specwf1()))
            print("Specw f2       = " +repr(get_c_specwf2()))
            print("Specw F3       = " +repr(get_c_specwf3()))
            print("OffF1          = " +repr(get_c_offsf1()))
            print("OffF2          = " +repr(get_c_offsf2()))
            print("OffF3          = " +repr(get_c_offsf3()))
            print("Max. Abs.      = " +repr(get_c_absmax()))


    def test_extb(self):
    # Test for peak extraction using backward linear prediction coeff 
    # calculated by LPSVD
        dim(1)
    #    unit( "h") 
        order( 10)
        simu( 1000, 512, 2, 1.0, 200, 1.0, 0., 1.0, 600, 1.0, 0., 0.0)
        dt2svd( 512)
        svd2ar( 2)
        ar2rt( 2)
        rtclean( 2)
        rtlist( 2,1,get_nrt())
        rtinv( 2)
        rtlist( 1,1,get_nrt())
        rt2pk( 512, 1, 0.)
        pklist(1,get_npk1d())
    #
    #Number of Root          2
    #backward roots of PE polynome
    #       -.809826
    #        .588373
    #        .309326
    #       -.952008
    #
    #Number of Root          2
    #forward roots of PE polynome
    #       -.808208
    #       -.587198
    #        .308708
    #        .950106
    #
    #Peaks computed with linear prediction
    #Number of 1D peaks :          2
    #Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase
    #     1   601.174    6.13590        1.00000          .318      .000
    #     2   200.391    6.13593       1.000000          .318      .000








    def test_extm(self):
    # Test for peak extraction by matching bacward and forward
    # roots
        dim(1)
        # unit( "h") 
        order( 10)
        simu( 1000, 512, 2, 1.0, 200, 1.0, 0., 1.0, 600, 1.0, 0., 0.0)
        dt2svd( 512)
        svd2ar( 1)
        ar2rt( 1)
        svd2ar( 2)
        ar2rt( 2)
        rtmatch( 3)
        rtlist( 3,1,get_nrt())
        rt2pk( 512, 2, 0.)
        pklist(1,get_npk1d())
    #
    #Number of Root          3
    #matched roots of PE polynome
    #        .308708
    #        .950106
    #       -.808208
    #       -.587198
    #        .990458
    #       -.161216
    #
    #Peaks computed with linear prediction
    #Number of 1D peaks :          3
    #Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase
    #     1   200.391    6.13592       1.000000          .318      .000
    #     2   601.174    6.13592        1.00000          .318      .000
    #     3   976.226      .000
    #       Linked to peak #     0







    """

    def fit():
    # this one test linefitting commands
        read( "fit_test")
        evaln( 33, 142); addbase( get_shidt()); com_max()
        zoom( 1, 636, 807)
        minimax(10*get_noise(), geta_max(2))
        peak()
        report_result( (get_npk1d == 21), "peak")

        zoom( 1, 643, 687)
        iter( 10)
        linefit( "lorentz")
        a = (get_chi2() < 37.4)
        linefit( "gauss")
        report_result( (a and (chi2 < 11.2)), "Linefit")

        zoom( 0)


        n = 100
        sp1 = 0;  sp2 = 0;  sdp1 = 0;  sdp2 = 0;  schi = 0
        print "Testing Chi2 law"
        refmacro 0
        for i = 1 to n
        #   print ("trial"#i)
           dim(1) chsize(50
           one() tm 50 50  put tab                # create a dump TAB[]
           one() mult 1.5 specw 1000 em 100       # create a dump data-set
           iter 10   miniter 10
           addnoise .02 123

           p1 := 1
           p2 := 0.1

          fitexp 1   # note that previous fit was x/p2 # while fitexp is x*p2
        #  print (p1#p2#dp1#dp2#chi2)
          sp1 = (%+p1)
          sdp1 = (%+dp1)
          sp2 = (%+p2)
          sdp2 = (%+dp2)
          schi = (%+chi2)
        endfor
        p1 = (sp1/n)
        p2 = (sp2/n)
        dp1 = (sdp1/n)
        dp2 = (sdp2/n)
        # if (GRAPH) .def macro/showexp 'p1*exp(-x*p2)'

        # print (p1#p2#dp1#dp2#schi/n)
        # print (p1/1.756#p2/7.824#dp1/0.03731#dp2/0.2170#schi/n)

        tol = (0.01)

        a = (math.fabs(1-p1/1.756)<tol)
        a = (a and (math.fabs(1-p2/7.824)<tol))
        a = (a and (math.fabs(1-dp1/0.03731)<tol))
        a = (a and (math.fabs(1-dp2/0.2170)<tol))
        a = (a and ((schi/n) < (55+1/sqrt(n)) ))
        report a "fitexp"

        p1 := -1
        p2 := 1
        MINIMIZE 'sin(p1)*sin(p2)' 2
        a = (math.fabs(p1+2*atan(1))<0.05)
        a = (a and (math.fabs(p2-2*atan(1))<0.05))
        a = (a and (math.fabs(chi2+1) < 1e-3))
        report a "minimize"

        unp1  unp2 undp1  undp2
        if (GRAPH) refmacro 1
    """
    def test_ft_base(self):
    # this one test basic reading and processing of 2Ds

        f = 'hsqc_grad'
        try:
            read(f)      # on sorbine
        except:
            print('** data/ directory is missing, test will not be executed **')

        com_max()
        self.report_result(get_dim() == 2 and get_si1_2d() == 256 and get_si2_2d() == 2048 and math.fabs(geta_max(1)-45705)<10.0, "Loading of data")

        com_sin(0.2, 'f12') 		# apodise
        ft_sim()    	# process in "f2" (simultaneous (complex) sampling here)
        chsize(2*get_si1_2d(), get_si2_2d())		# zerofill in F1
        ft_sh_tppi()	# process in "f1" (States-Haberkorn + TPPI here)
        phase(119.2, -80.0, 'f2')

        com_max()
        self.report_result((math.fabs(geta_max(1)-40154644.)<40.0), "Basic Fourier Transform and phasing")

        real('f12')
        com_max()
        report_result(math.fabs(geta_max(1)-39193776.)<40.0, "Real")

        row(86)   				# select a row
        dim(1); extract(584, 922); com_max()
        a = (math.fabs(geta_max(1)-4283381.)<10.0)
        dim(2); col( 819); dim(1)
        com_max()
        report  (math.fabs(geta_max(1)-7945463.)<8.0 and a, "dim , row, col, extract")

        dim(2); proj("f1", "s"); dim(1)
        evaln(33, 142); com_max()
        report (math.fabs(geta_max(2)-361132.4)<4.0 and math.fabs(noise-81255)<1.0, "proj evaln")

        writec( "fit_test")


    def test_inout(self):
    # to test for reading and writing files.
    # not in the standard test because the file used are not distributed
    #
    # standard file are not tested, because they are in cache

        dim(1)
        readv("data_test/hpf576hpresat2611af.fid/fid")
        com_max()
        self.report_result( ((geta_max(1) == 219632) and (get_si1_1d() == 6848)), 'ReadV 1D')

        dim(2)
        readv("data_test/hpf576hnoe1502611af.fid/fid")
        com_max()
        self.report_result( ((geta_max(1) == 111982) and (get_si1_2d() == 768) and (si2_2d == 3008)), 'ReadV 2D')

    """
    The following tests are are finished yet !


    def linpred
    # Test for forward linear prediction calculated by LPSVD
    # T.Malliavin first version
    # G.Salnikov jan-1996 : added peak sorting

        dim(1)
        order( 10)
        simu 6024 512 2 1.0 3123 1.0 0. 1.0 5001 1.0 0. 0.0
        dt->svd 512
    svd->ar 1
    ar->dt 1k 1
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - .95968282)<0.001)
    b = (math.fabs(shift - .226734569878E-03)<0.001)
    report (a and b) "Forward Linear Prediction by LPSVD"
    # forward linear prediction coeff calculated by LPSVD
    #Number of AR coefficients :         10
    #Forward coeff. of PE polynome
    #        .032423
    #        .218044
    #       -.135395
    #       -.069564
    #       -.022323
    #        .031589
    #        .083497
    #       -.560398
    #        .330764
    #        .119651
    #       -.255565
    #       -.082291
    #        .069599
    #        .047879
    #        .009635
    #        .261716
    #        .309227
    #       -.071307
    #       -.189211
    #        .000804

    # Test for backward linear prediction by LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 2
    ar->dt 1k 2
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - 1.00010108948)<0.001)
    b = (math.fabs(shift - .145896791946E-02)<0.001)
    report (a and b) "Backward Linear Prediction by LPSVD"

    # backward linear prediction coeff calculated by LPSVD
    #Number of AR coefficients :         10
    #Backward coeff. of PE polynome
    #        .211119
    #       -.018458
    #       -.129454
    #        .013135
    #        .123079
    #       -.092347
    #        .236444
    #        .761579
    #        .149035
    #        .181000
    #       -.181964
    #        .033043
    #        .164967
    #        .020283
    #       -.073370
    #       -.050815
    #        .152297
    #        .031029
    #       -.123679
    #       -.005283

    # Test for forward linear prediction and spectrum calculation by the Burg method
    dim(1)
    order 10
    simu 6024 512 1 1.0 3123 1.0 0. 0.0
    dt->ar
    ar->dt 1k 1
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - .684671044350)<0.001)
    b = (math.fabs(shift - .664793420583E-03)<0.001)
    ar->sp
    com_max()
    c = (math.fabs(geta_max(1) - .312438824039E-04)<0.001)
    d = (math.fabs(geta_max(2) - .178954522312E-07)<0.001)
    report (a and b and c and d) "Forward Linear Prediction and spectrum calculation by the Burg method"
    # Linear prediction coefficients calculated by the Burg method
    #Number of AR coefficients :         10
    #Forward coeff. of PE polynome
    #       1.298755
    #        .151038
    #        .205728
    #        .048508
    #       -.028264
    #       -.010233
    #        .012509
    #        .006233
    #       -.006552
    #       -.004286
    #        .003542
    #        .002954
    #       -.001666
    #       -.001754
    #        .000144
    #        .000191
    #        .004596
    #        .007862
    #        .018509
    #        .042231

    # Test for spectral analysis by forward LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 1
    ar->rt 1 
    rt->pk 512 1 0.0
    unit h
    pkclean Below 0.001
    for i = 1 to npk1d
    f[i] := pk1d_f[i]
    w[i] := pk1d_w[i]
    endfor
    for i = 1 to (npk1d-1)
      if (f[i] > f[i+1]) then
        for j = i to 1 step -1
          if (f[j] > f[j+1]) then
    	 k = f[j+1]
    	 f[j+1] = f[j]
    	 f[j] = k
    	 k = w[j+1]
    	 w[j+1] = w[j]
    	 w[j] = k
          endif
        endfor
      endif
    endfor
    a = (math.fabs(f[1] - 86.9482116699) < 0.001)
    b = (math.fabs(f[2] - 246.565734863) < 0.001)
    c = (math.fabs(w[1] - .270543918014E-01) < 0.001)
    d = (math.fabs(w[2] - .270575433969E-01) < 0.001)
    report (a and b and c and d) "Spectral analysis by Forward LPSVD"
    #Peaks computed with linear prediction
    #Number of 1D peaks :         10
    #Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase  Type
    #     1  4273.871    .112823E-06    .904885E-06    94.365    58.986 unknown
    #     2  1955.547    .479791E-07    .853707E-06   209.350   152.705 unknown
    #     3   447.116    .327822E-07    .242444E-06    87.014    98.625 unknown
    #     4  2115.433    .491742E-07    .184734E-05   442.003   -58.959 unknown
    #     5  1101.549    .631592E-08    .297399E-06   554.010     6.059 unknown
    #     6  5697.313    .139044E-06    .734033E-05   621.123   124.380 unknown
    #     7  5010.787    36.9626        1.00000          .318      .000 unknown
    #     8  3129.112    36.9584        1.00000          .318      .000 unknown
    #     9  3575.401    .111505E-06    .785088E-06    82.840    41.604 unknown
    #    10  5886.629    .101792E-06    .952602E-05  1101.066   -90.490 unknown
    pkselect 1 2 0
    for i = 1 to npk1d
    f[i] := pk1d_f[i]
    w[i] := pk1d_w[i]
    endfor
    for i = 1 to (npk1d-1)
      if (f[i] > f[i+1]) then
        for j = i to 1 step -1
          if (f[j] > f[j+1]) then
    	 k = f[j+1]
    	 f[j+1] = f[j]
    	 f[j] = k
    	 k = w[j+1]
    	 w[j+1] = w[j]
    	 w[j] = k
          endif
        endfor
      endif
    endfor
    a = (math.fabs(f[1] - 86.9482116699) < 0.5)
    b = (math.fabs(w[2] - 0.270575117320E-01) < 0.5)
    report (a and b) "PKSELECT"
    #Peaks computed with linear prediction
    #Number of 1D peaks :          2
    # Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase  Type
    #     1  1955.547    .479791E-07    .853707E-06   209.350   152.705 Gauss
    #     2  2115.433    .491742E-07    .184734E-05   442.003   -58.959 Gauss


    # Test for spectral analysis by backward LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 2
    ar->rt 2
    rtclean 2
    rtinv 2 
    rt->pk 512 1 0.0
    unit h
    for i = 1 to npk1d
    f[i] := pk1d_f[i]
    w[i] := pk1d_w[i]
    endfor
    for i = 1 to (npk1d-1)
      if (f[i] > f[i+1]) then
        for j = i to 1 step -1
          if (f[j] > f[j+1]) then
    	 k = f[j+1]
    	 f[j+1] = f[j]
    	 f[j] = k
    	 k = w[j+1]
    	 w[j+1] = w[j]
    	 w[j] = k
          endif
        endfor
      endif
    endfor
    a = (math.fabs(f[2] - 246.565734863) < 0.001)
    b = (math.fabs(f[1] - 86.9482116699) < 0.001)
    c = (math.fabs(w[2] - .270543918014E-01) < 0.001)
    d = (math.fabs(w[1] - .270575433969E-01) < 0.001)
    e = (nrt == 2)
    report (a and b and c and d &e) "Spectral analysis by Backward LPSVD"
    #Number of 1D peaks :          2
    #Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase  Type
    #     1  3129.112    36.9573        1.00000          .318      .000 unknown
    #     2  5010.787    36.9623        1.00000          .318      .000 unknown

    # Test for phase-minimum bacward and forward linear prediction by LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 2
    ar->rt 2
    rtclean 2
    rtinv 2 
    rt->ar 3
    #Number of AR coefficients :          2
    #Forward coeff. of PE polynome
    #        .510480
    #        .991116
    #       -.580482
    #        .813865
    #Backward coeff. of PE polynome
    #        .510650
    #       -.991445
    #       -.580867
    #       -.814406
    ar->dt 1k 1
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - .959682941437)<0.001)
    b = (math.fabs(shift - .226510470384E-03)<0.001)
    report (a and b) "Phase-minimum forward Linear Prediction by LPSVD"
    chsize(512
    ar->dt 1k 2
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - 1.00010097027)<0.001)
    b = (math.fabs(shift - .145897129551E-02)<0.001)
    report (a and b) "Phase-minimum backward Linear Prediction by LPSVD"

    # Test for data recalculation from SVD decomposition
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    mult -1 put data mult -1
    dt->svd 512
    #number of SVD :         10
    #Singular values
    #      49.514317
    #        .000000
    #      47.397547
    #        .000000
    #        .000176
    #        .000000
    #        .000167
    #        .000000
    #        .000137
    #        .000000
    #        .000080
    #        .000000
    #        .000064
    #        .000000
    #        .000056
    #        .000000
    #        .000054
    #        .000000
    #        .000046
    #        .000000
    svdclean2 5.0 y
    c = (nsvd == 2)
    svdclean1 2 y
    svd->dt
    adddata
    com_max()
    a = (math.fabs(geta_max(1) - .160932540894E-04) < 0.001)
    b = (math.fabs(geta_max(2) + .217556953430E-04) < 0.001)
    report (a and b &c) "Data back-calculation from SVD decomposition"

    # Test for forward linear prediction calculated by LPSVD using roots reflected into the unit-circle
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 1
    ar->rt 1
    rtreflect 1
    rt->ar 1
    ar->dt 1k 1
    itype 0 evaln 1 1k itype 1
    a = (math.fabs(noise - .95968062)<0.001)
    b = (math.fabs(shift - .22652098E-03)<0.001)
    report (a and b) "Forward Linear Prediction by LPSVD using reflected roots"
    # forward linear prediction coeff calculated by LPSVD
    #Number of AR coefficients :         10
    #Forward coeff. of PE polynome
    #        .032423
    #        .218044
    #       -.135395
    #       -.069564
    #       -.022323
    #        .031589
    #        .083497
    #       -.560398
    #        .330764
    #        .119651
    #       -.255565
    #       -.082291
    #        .069599
    #        .047879
    #        .009635
    #        .261716
    #        .309227
    #       -.071307
    #       -.189211
    #        .000804

    # Test for spectral analysis by backward LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 3
    ar->rt 3
    rtmatch 2 
    rt->pk 512 2 0.0
    unit h
    for i = 1 to npk1d
    f[i] := pk1d_f[i]
    w[i] := pk1d_w[i]
    endfor
    for i = 1 to (npk1d-1)
      if (f[i] > f[i+1]) then
        for j = i to 1 step -1
          if (f[j] > f[j+1]) then
    	 k = f[j+1]
    	 f[j+1] = f[j]
    	 f[j] = k
    	 k = w[j+1]
    	 w[j+1] = w[j]
    	 w[j] = k
          endif
        endfor
      endif
    endfor
    a = (math.fabs(f[2] - 246.565734863) < 0.001)
    b = (math.fabs(f[1] - 86.9482116699) < 0.001)
    c = (math.fabs(w[2] - .270579215139E-01) < 0.001)
    d = (math.fabs(w[1] - .270544942468E-01) < 0.001)
    e = (nrt == 2)
    report (a and b and c and d &e) "Spectral analysis by Matched LPSVD"
    #Peaks computed with linear prediction
    #Number of 1D peaks :          2
    #Peak coordinates in HZ       widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase  Type
    #     1  3129.112    36.9579        1.00000          .318      .000 Gauss
    #     2  5010.787    36.9624        1.00000          .318      .000 Gauss

    # Test for spectral analysis by forward LPSVD
    dim(1)
    order 10
    simu 6024 512 2 1.0 3123 1.0 0. \
    1.0 5001 1.0 0. 0.0
    dt->svd 512
    svd->ar 1
    ar->rt 1
    rtselect 1 1 2 3 4 5 6 7 8 9 10 0
    a = (nrt == 10)
    rt->pk 512 1 0.
    pkclean Below 0.001
    for i = 1 to npk1d
    f[i] := pk1d_f[i]
    w[i] := pk1d_w[i]
    endfor
    for i = 1 to (npk1d-1)
      if (f[i] > f[i+1]) then
        for j = i to 1 step -1
          if (f[j] > f[j+1]) then
    	 k = f[j+1]
    	 f[j+1] = f[j]
    	 f[j] = k
    	 k = w[j+1]
    	 w[j+1] = w[j]
    	 w[j] = k
          endif
        endfor
      endif
    endfor
    b = (math.fabs(f[1] - 86.9482116699) < 0.5)
    c = (math.fabs(f[2] - 246.565734863) < 0.5)
    d = (math.fabs(w[1] - 0.270543936640E-01) < 0.5)
    e = (math.fabs(w[2] - 0.270575117320E-01) < 0.5)
    report (a and b and c and d &e) "RTSELECT"
    #Number of Root          3
    #forward roots of PE polynome
    #       -.247805
    #       -.919160
    #        .853978
    #        .428907
    #        .678869
    #       -.249615
    #Peaks computed with linear prediction
    #Number of 1D peaks :          3
    #Peak coordinates in INDEX    widthes in Hz
    #   Index   Coord          Ampli         Volume       Width     Phase  Type
    #     1   149.459    .659790E-02    .529180E-01    94.365    65.051 unknown
    #     2   474.072    .150557E-01    .111346        87.014    46.328 unknown
    #     3    28.712    .191453E-01    1.01070       621.123   -41.029 unknown

    # Test for forward linear prediction and spectrum calculation by the Burg method (BURG command)
    dim(1)
    order 10
    simu 6024 512 1 1.0 3123 1.0 0. 0.0
     if (config_os s= "SOLARIS") then
      print "SORRY, BURG is not implemented on Solaris yet"
     else
      burg 1k
      itype 0 evaln 1 1k itype 1
      a = (math.fabs(noise - .684670567513)<0.001)
      b = (math.fabs(shift - .664792780299E-03)<0.001)
      report (a and b) "BURG command"
     endif

    # Linear prediction coefficients calculated by the Burg method
    #Number of AR coefficients :         10
    #Forward coeff. of PE polynome
    #       1.298755
    #        .151038
    #        .205728
    #        .048508
    #       -.028264
    #       -.010233
    #        .012509
    #        .006233
    #       -.006552
    #       -.004286
    #        .003542
    #        .002954
    #       -.001666
    #       -.001754
    #        .000144
    #        .000191
    #        .004596
    #        .007862
    #        .018509
    #        .042231
    """
    def test_maxentropy(self):
        """
    # This test should be considerably extended.
    # for the moment, only Inverse Laplace/Tabulated is tested
    # the minimum would be to add 2D FT Maxent
        """
        dim(1); itype( 0)
        dmin( .001)
        dmax( 1)
        # create dump tab[]
        chsize(30)
        for i in range(30):
          setval( i+1, 0.1*math.exp(i+1/3))
        put("tab")

        ntrial=4
        seed=127
        for i in range(ntrial):  # try several times, as it may fail (because of addnoise)
            #create asym damping distrib
            chsize(15); one(); mult( 1000); tm( 1, 1)
            chsize(30); reverse(); chsize(40)
            tlaplace()
            ns=0.001*val1d(1)
            addnoise (ns, seed)
            seed=seed+11
            noise(ns)

            #gifa
            algo(11)
            iter( 200); miniter( 4); ndisp( 50);  lambsp( 4)
            invtlap( 100)
            iter( 1000)
            invtlapcont()

            # I get on my ref machine (linux RH 6) :
            # noise = 115 - 118  shift = 80.1 - 80.2
            # #68 = 356 - 376
            # chi2 = 1.2 - 2.4

            a = (get_iterdone() > 10 and get_iterdone() <= 1000 and get_chi2() < 3)  # iterdone() == 10 seems unlikely !
            evaln( 1, 100)
            a = (a and (get_shift()<80.5) and (get_noise()<125) and val1d(68)>330)
            writec("test-r-%d.gs1"%(i))
        self.report_result( a, 'Inverse Laplace Transform (tried %i time)'%(ntrial))

        if ( not(a and get_chi2()  < 5)):          # if chi2 is off but not so badly
                    print("iterdone =",get_iterdone(),"    CHI2 =",get_chi2(), "    SHIFT=",get_shift(), "     NOISE=",get_noise(),"   VAL(68)=",val1d(68))
                    print('the Inverse Laplace Transform fails because of a convergence slower than normal')
                    print('You can still use this module, but will have not optimal convergence speed')

    def test_memory(self):
        """check for dynamic memory size
        """
    #    l = (1,4,16,64,256,1024,4096,16*1024,64*1024)
        l = (1,4,16,64,256,1024,4096,16*1024)
        for i in l:
            dim(1)
            n=i*1024
            chsize(n)
            self.report_result(1,"1D buffer size: "+str(n)) 
        dim(1)
        chsize(1024)
        for i in l:
            dim(2)
            chsize(1024,i)
            self.report_result(1,"2D buffer size: 1024 x "+str(i)) 
        dim(2)
        chsize(128,128)
        for i in l:
            dim(3)
            chsize(32,32,i)
            self.report_result(1,"3D buffer size: 32 x 32 x "+str(i)) 
        dim(3)
        chsize(16,16,16)

    def test_pkhandle(self):
        """
    # testing PKSYM in 2D, unit Hertz and PPM (using SIMUN)
        """
        # create dummy spec and peak pick
        dim(2); specw( 5283.0, 4230.); freq( 600.13, 600.13, 600.); offset(500., 500.)
        chsize(512, 1024); zero(); unit( "h")
        simun( 1, 1.0, 3461, 999, 2.1, 1.5, 0., 0.)
        simun( 1, 1.0, 5000, 601, 2.1, 1.5, 0., 0.)
        sin( 0, 'f12'); ft( 'f12'); com_max()
        minimax( 5000, geta_max(1)*2)
        real( 'f12')
        zoom( 0)
        peak( 0)

        if (get_npk2d() == 2):
          a = (math.fabs(geta_pk2d_f1f(1) - 114)<0.01)
          b = (math.fabs(geta_pk2d_f1f(2) - 39)<0.01)
          c = (math.fabs(geta_pk2d_f2f(1) - 453)<0.01)
          d = (math.fabs(geta_pk2d_f2f(2) - 501)<0.01)
          e = (math.fabs(geta_pk2d_a(1) - 9228.83)<1.0)
          tt = (a and b and c and d and e)
        else:
          tt = (1 == 0)
        self.report_result( tt, "PEAK in 2D")

        pksym( 1, 1.0)
        a = (math.fabs(geta_pk2d_f1f(3) - 232.98 )<0.01)
        b = (math.fabs(geta_pk2d_f2f(3) - 155.804)<0.01)
        c = (math.fabs(geta_pk2d_f1f(4) - 252.196)<0.01)
        d = (math.fabs(geta_pk2d_f2f(4) + 31.536 )<0.01)
        self.report_result( (a and b and c and d), "PKSYM 2D with option add")

        pkselect( 1, 2, 3, 0)
        pksym( 0, 1.0)
        a = (npk2d == 2)
        b = (math.fabs(geta_pk2d_f1f(1) - 114)<0.01)
        c = (math.fabs(geta_pk2d_f1f(2) - 232.980)<0.01)
        self.report_result( (a and b and c), "PKSELECT - PKSYM 2D with option remove")

    """
        pkhandle_test_file = (gifa_tmpdir//file_separator//'pkhandle_test_file')
        tab2pk pkhandle_test_file new
        pkclear
        report (npk2d == 0)  "PKCLEAR"

        pk2tab pkhandle_test_file new
        a = (npk2d == 2)
        pkclear
        transpose peak 0
        pk2tab pkhandle_test_file append
        b = (math.fabs(pk2d_f1f[1] - 501)<0.01)
        c = (math.fabs(pk2d_f1f[3] - 453.03)<0.01)
        report (a and npk2d == 4 and b and c ) "tab2pk and pk2tab"

        rm (pkhandle_test_file // dbext)

        # testing PKPROJ in 2D, unit Hertz (using SIMUN)

        dim(2) specw 5283.0 4230. freq 600.13 600.13 600. off500. 500.
        chsize(512 1k zero() unit h

        simun 1 1.0 3461 999 2.1 1.5 0. 0.
        simun 1 1.0 5000 601 2.1 1.5 0. 0.
        sin 0. f12 ft f12 com_max()
        minicom_max() 5000 (geta_max(1)*2)
        real f12 peak 0
        dim(1) specw 5283.0 freq 600.13 600. off500. chsize(256 dim(2)
        pkproj "f2" s 1.0
        c = (math.fabs(pk1d_f[1] - 114)<0.01)
        d = (math.fabs(pk1d_f[2] - 39)<0.01)
        report (c and d) "PKPROJ 2D along F2"
        dim(1) specw 4230.0 freq 600.13 600. off500. chsize(512 dim(2)
        pkproj "f1" s 1.0
        a = (math.fabs(pk1d_f[1] - 453)<0.01)
        b = (math.fabs(pk1d_f[2] - 501)<0.01)
        report (a and b) "PKPROJ 2D along F1"
        dim(1) pkclear
        dim(2) pkclear
    """




    """
    def simu_t
    # testing SIMUN function in 1D
    dim(1) specw 6283.185 freq 600.13 600. off0.
    chsize(1k zero() unit i
    simun 1.0 221 1.0 0.0 
    ft com_max()
    a = (math.fabs(geta_max(1) - 419.92673)<0.001)
    b = (math.fabs(geta_max(2) + 140.58662)<0.001)
    minicom_max() (geta_max(1)/2) (geta_max(1)*2)
    real peak
    c = (math.fabs(pk1d_f[1] - 111.0)<0.001)
    report (a and b and c) "SIMUN 1D with index unit"

    chsize(1k zero() unit h
    simun 1.0 3461 1.0 0.0 
    simun 1.0 2161 1.0 0.0 
    sin 0. ft com_max()
    a = (math.fabs(geta_max(1) - 296.48993)<0.001)
    b = (math.fabs(geta_max(2) + 109.73332)<0.001)
    minicom_max() (geta_max(1)/2) (geta_max(1)*2)
    real peak
    c = (math.fabs(pk1d_f[1] - 231.0)<0.001)
    d = (math.fabs(pk1d_f[2] - 337.0)<0.001)
    report (a and b and c and d) "SIMUN 1D with hertz unit"

    dim(1) specw 6283.185 freq 600.13 600.
    chsize(1k zero() unit p off500.
    simun 1.0 7.0 20.0 0.0 
    simun 1.0 3.0 33.0 0.0 
    ft com_max()
    a = (math.fabs(geta_max(1) - 73.347198)<0.001)
    b = (math.fabs(geta_max(2) + 44.897629)<0.001)
    minicom_max() (geta_max(1)/2) (geta_max(1)*2) real peak
    evaln 10 50 linefit %
    c = (math.fabs(pk1d_f[1] - 211.496780396)<0.001)
    d = (math.fabs(pk1d_f[2] - 407.066345215)<0.001)
    e = (math.fabs(pk1d_w[1] - 1.63042700291)<0.001)
    f = (math.fabs(pk1d_w[2] - 2.68996524811)<0.001)
    report (a and b and c and d and e and f) "SIMUN 1D with PPM unit, off= 500 and linefit"

    # testing SIMUN function in 2D, unit index, amplitude modulation
    dim(2) specw 7812.5 1893.939 freq 600.13 150.9 600. off500. 500.
    chsize(512 1k zero() unit i
    simun 1 1.0 241 670 2.1 1.5 0. 0.
    simun 1 1.0 417 365 2.1 1.5 0. 0. 
    sin 0. f12 ft f12 com_max()
    a = (math.fabs(geta_max(1) - 24396.150)<0.1)
    b = (math.fabs(geta_max(2) + 18912.082)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real f12 peak 0
    c = (math.fabs(pk2d_f2f[1] - 183.0)<0.1)
    d = (math.fabs(pk2d_f1f[1] - 209.0)<0.1)
    e = (math.fabs(pk2d_f2f[2] - 336.0)<0.1)
    f = (math.fabs(pk2d_f1f[2] - 121.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with index unit and ampl. modul."

    # testing SIMUN function in 2D, unit Hertz,amplitude modulation
    dim(2) specw 7812.5 1893.939 freq 600.13 150.9 600. off500. 500.
    chsize(512 1k zero() unit h
    simun 1 1.0 3461 999 2.1 1.5 0. 0. 
    simun 1 1.0 7023 601 2.1 1.5 0. 0. 
    #imun 1 1.0 6001 999 2.1 1.5 0. 0. 
    #imun 1 1.0 2541 227 2.1 1.5 0. 0. 
    sin 0. f12 ft f12 com_max()
    a = (math.fabs(geta_max(1) - 31308.023)<0.1)
    b = (math.fabs(geta_max(2) + 14715.346)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real f12 peak 0
    c = (math.fabs(pk2d_f2f[1] - 378.0)<0.1)
    d = (math.fabs(pk2d_f1f[1] - 160.0)<0.1)
    e = (math.fabs(pk2d_f2f[2] - 486.0)<0.1)
    f = (math.fabs(pk2d_f1f[2] - 43.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with Hertz unit and ampl. modul."

    # testing SIMUN function in 2D, unit PPM, amplitude modulation
    chsize(512 1k zero() unit p
    simun 1 1.0 36 2.5 2.1 1.5 0. 0. 
    simun 1 1.0 43 1.2 2.1 1.5 0. 0. 
    sin 0. f12 ft f12 com_max()
    a = (math.fabs(geta_max(1) - 15827.653)<0.1)
    b = (math.fabs(geta_max(2) + 16102.342)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real f12 peak 0
    c = (math.fabs(pk2d_f2f[1] - 243.0)<0.1)
    d = (math.fabs(pk2d_f1f[1] - 95.0)<0.1)
    e = (math.fabs(pk2d_f2f[2] - 454.0)<0.1)
    f = (math.fabs(pk2d_f1f[2] - 61.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with PPM unit and ampl. modul."

    # testing SIMUN function in 2D, unit index, phase modulation
    dim(2) specw 7812.5 1893.939 freq 600.13 150.9 600. off300. 200.
    chsize(512 1k zero() unit i
    simun 0 1.0 241 670 2.1 1.5 0. 0. 
    simun 0 1.0 217 365 2.1 1.5 0. 0. 
    sin 0.5 f12 ft "f2" flip ft "f1" flop com_max()
    a = (math.fabs(geta_max(1) - 35222.922)<0.1)
    b = (math.fabs(geta_max(2) + 41904.309)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real "f2" peak 0
    c = (math.fabs(pk2d_f2f[2] - 184.0)<0.1)
    d = (math.fabs(pk2d_f1f[2] - 217.0)<0.1)
    e = (math.fabs(pk2d_f2f[4] - 336.0)<0.1)
    f = (math.fabs(pk2d_f1f[4] - 241.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with index unit and phase modul."

    # testing SIMUN function in 2D, unit Hertz, phase modulation
    dim(2) specw 7812.5 1893.939 freq 600.13 150.9 600. off300. 200.
    chsize(512 1k zero() unit h
    simun 0 1.0 6001 999 2.1 1.5 0. 0. 
    simun 0 1.0 2541 227 2.1 1.5 0. 0. 
    sin 0.5 f12 ft "f2" flip ft "f1" flop com_max()
    a = (math.fabs(geta_max(1) - 43453.328)<0.1)
    b = (math.fabs(geta_max(2) + 42027.289)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real "f2" peak 0
    c = (math.fabs(pk2d_f2f[2] - 297.0)<0.1)
    d = (math.fabs(pk2d_f1f[2] - 139.0)<0.1)
    e = (math.fabs(pk2d_f2f[4] - 506.0)<0.1)
    f = (math.fabs(pk2d_f1f[4] - 366.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with Hertz unit and phase modul."

    # testing SIMUN function in 2D, unit PPM, phase modulation
    dim(2) specw 7812.5 1893.939 freq 600.13 150.9 600. off530. 101.
    chsize(512 1k zero() unit p
    simun 0 1.0 36 2.5 2.1 1.5 0. 0. 
    simun 0 1.0 43 1.2 2.1 1.5 0. 0. 
    sin 0.5 f12 ft "f2" flip ft "f1" flop com_max()
    a = (math.fabs(geta_max(1) - 46297.816)<0.1)
    b = (math.fabs(geta_max(2) +30157.898)<0.1)
    minicom_max() 5000 (geta_max(1)*2)
    real "f2" peak 0
    c = (math.fabs(pk2d_f2f[1] - 134.0)<0.1)
    d = (math.fabs(pk2d_f1f[1] - 192.0)<0.1)
    e = (math.fabs(pk2d_f2f[2] - 345.0)<0.1)
    f = (math.fabs(pk2d_f1f[2] - 123.0)<0.1)
    report (a and b and c and d and e and f) "SIMUN 2D with PPM unit and phase modul."

    # testing SIMUN function in 3D, unit index
    dim 3 specw 100.0 100.0 100.0 freq 400.0 400.0 400.0 400. off8. 0. 5.
    chsize(64 32 128 zero() unit i
    simun 1.0 21. 14. 56. 0.2 0.2 0.2 0. 0. 0. 
    simun 1.0 10. 27. 99. 0.2 0.2 0.2 0. 0. 0. 
    sin 0. f123 ft f123 com_max()
    a = (math.fabs(geta_max(1) - 3064.7778)<0.01)
    b = (math.fabs(geta_max(2) + 2398.55639648)<0.01)
    minicom_max() 100 3100
    real f123 peak 0
    c = (math.fabs(pk3d_f3f[1] - 50.0)<0.01)
    d = (math.fabs(pk3d_f2f[1] - 14.0)<0.01)
    e = (math.fabs(pk3d_f1f[1] - 6.0)<0.01)
    f = (math.fabs(pk3d_f3f[2] - 29.0)<0.01)
    g = (math.fabs(pk3d_f2f[2] - 8.0)<0.01)
    h = (math.fabs(pk3d_f1f[2] - 11.0)<0.01)
    report (a and b and c and d and e and f and g and h) "SIMUN 3D with index unit"

    # testing SIMUN function in 3D, unit Hertz
    dim 3 specw 100.0 100.0 100.0 freq 400.0 400.0 400.0 400. off8. 0. 5.
    chsize(64 32 128 zero() unit h
    simun 1.0 50. 29. 56. 0.2 0.2 0.2 0. 0. 0. 
    simun 1.0 31. 73. 82. 0.2 0.2 0.2 0. 0. 0. 
    sin 0. f123 ft f123 com_max()
    a = (math.fabs(geta_max(1) - 1940.9633)<0.01)
    b = (math.fabs(geta_max(2) + 1932.8833)<0.01)
    minicom_max() 1000 3000
    real f123 peak 0
    c = (math.fabs(pk3d_f3f[1] - 32.0)<0.01)
    d = (math.fabs(pk3d_f2f[1] - 12.0)<0.01)
    e = (math.fabs(pk3d_f1f[1] - 20.0)<0.01)
    f = (math.fabs(pk3d_f3f[2] - 16.0)<0.01)
    g = (math.fabs(pk3d_f2f[2] - 5.0)<0.01)
    h = (math.fabs(pk3d_f1f[2] - 26.0)<0.01)
    report (a and b and c and d and e and f and g and h) "SIMUN 3D with Hertz unit"

    # testing SIMUN function in 3D, unit PPM
    dim 3 specw 100.0 100.0 100.0 freq 400.0 400.0 400.0 400. off8. 0. 5.
    chsize(64 32 128 zero() unit p
    simun 1.0 0.1 0.2 0.15 0.2 0.2 0.2 0. 0. 0. 
    simun 1.0 0.05 0.17 0.21 0.2 0.2 0.2 0. 0. 0. 
    sin 0. f123 ft f123 com_max()
    a = (math.fabs(geta_max(1) - 3703.9995)<0.01)
    b = (math.fabs(geta_max(2) + 3339.1870)<0.01)
    minicom_max() 1000 4000
    real f123 peak 0
    c = (math.fabs(pk3d_f3f[1] - 30.0)<0.01)
    d = (math.fabs(pk3d_f2f[1] - 4.0)<0.01)
    e = (math.fabs(pk3d_f1f[1] - 23.0)<0.01)
    f = (math.fabs(pk3d_f3f[2] - 14.0)<0.01)
    g = (math.fabs(pk3d_f2f[2] - 6.0)<0.01)
    h = (math.fabs(pk3d_f1f[2] - 29.0)<0.01)
    report (a and b and c and d and e and f and g and h) "SIMUN 3D with PPM unit"

    quit
    """

    def report_result(self,t,msg):
        """report test message

        used by test suite
        """
        l = len(msg)
        l = (max(1,80-l))
        if t :
                print("--- " , msg  , "."*(l-4) , " Ok")
        else: 
                print("--- " , msg  , "."*(l-8) , " FAILED")

    # def test_report(self,msg):
    #     """report test message
    # 
    #     used by test suite
    #     """
    #     l = len(msg)
    #     l = (max(1,80-l))
    #     print "\n### " , msg  , "#"*l

    # def main():
    #     """docstring for dotest"""
    # 
    #     report("Basic")
    #     basic()
    #     # report("pkhandle")
    #     # pkhandle()
    #     report("Benchmark")
    #     benchmark()
    #     #report("Maxentropy")
    #     #maxentropy()
    #     report("Memory")
    #     memory()
    #     report("cache")
    #     cache()

    def new_func(x):
        print(x)

if __name__ == "__main__":
    unittest.main()

