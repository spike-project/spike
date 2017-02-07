#!/usr/bin/env python 
# encoding: utf-8

"""
    This file contains the definitions of all the NMR-active spins
    
    every thing is held into a big table:

 table[isotope-name]= spin,  naturalAbondance, magneticMoment, magnetogyricRatio, freq-at-100MHz,  quadrupoleMoment

 spin               in number of 1/2 (i.e spin == 1 means Spin is 1/2)
 naturalAbondance   in %
 magnetogyricRatio  in 10^7 rad/Tesla
 freq-at-100MHz     in MHz
 quadrupoleMoment   in fm2


Values (tentatively) from
ROBIN K. HARRIS, EDWIN D. BECKER, SONIA M. CABRAL DE MENEZES, ROBIN GOODFELLOW, AND PIERRE GRANGER
"NMR NOMENCLATURE.NUCLEAR SPIN PROPERTIES AND CONVENTIONS FOR CHEMICAL SHIFTS (IUPAC Recommendations 2001)"
Pure Appl.Chem., Vol.73, No.11, pp.1795-1818, 2001.
"""

from __future__ import print_function, division

__author__ = "Marc A. Delsuc <delsuc@igbmc.fr>"
#__date__ = "Oct 2009"  # MAD added  table4 
__date__ = "Jul 2010"   # MAD minor corrections using pylint



# table[isotope-name]= (spin,  naturalAbondance, magneticMoment, magnetogyricRatio, freq,  quadrupoleMoment)
table = {}
# dipolar
table["1H"] = (1, 99.9885, 4.837353570, 26.7522128, 100.000000, 0.0)
table["3H"] = (1, 0, 5.159714367, 28.5349779, 106.663974, 0.0)
table["3He"] = (1, 0.000137, -3.685154336, -20.3801587, 76.179437, 0.0)
table["13C"] = (1, 1.07, 1.216613, 6.728284, 25.145020, 0.0)
table["15N"] = (1, 0.368, -0.49049746, -2.71261804, 10.136767, 0.0)
table["19F"] = (1, 100, 4.553333, 25.18148, 94.094011, 0.0)
table["29Si"] = (1, 4.6832, -0.96179, -5.3190, 19.867187, 0.0)
table["31P"] = (1, 100, 1.95999, 10.8394, 40.480742, 0.0)
table["57Fe"] = (1, 2.119, 0.1569636, 0.8680624, 3.237778, 0.0)
table["77Se"] = (1, 7.63, 0.92677577, 5.1253857, 19.071513, 0.0)
table["89Y"] = (1, 100, -0.23801049, -1.3162791, 4.900198, 0.0)
table["103Rh"] = (1, 100, -0.1531, -0.8468, 3.186447, 0.0)
table["107Ag"] = (1, 51.839, -0.19689893, -1.0889181, 4.047819, 0.0)
table["109Ag"] = (1, 48.161, -0.22636279, -1.2518634, 4.653533, 0.0)
table["111Cd"] = (1, 12.80, -1.0303729, -5.6983131, 21.215480, 0.0)
table["113Cd"] = (1, 12.22, -1.0778568, -5.9609155, 22.193175, 0.0)
table["115Sn"] = (1, 0.34, -1.5915, -8.8013, 32.718749, 0.0)
table["117Sn"] = (1, 7.68, -1.73385, -9.58879, 35.632259, 0.0)
table["119Sn"] = (1, 8.59, -1.81394, -10.0317, 37.290632, 0.0)
table["123Te"] = (1, 0.89, -1.276431, -7.059098, 26.169742, 0.0)
table["125Te"] = (1, 7.07, -1.5389360, -8.5108404, 31.549769, 0.0)
table["129Xe"] = (1, 26.44, -1.347494, -7.452103, 27.810186, 0.0)
table["183W"] = (1, 14.31, 0.20400919, 1.1282403, 4.166387, 0.0)
table["187Os"] = (1, 1.96, 0.1119804, 0.6192895, 2.282331, 0.0)
table["195Pt"] = (1, 33.832, 1.0557, 5.8385, 21.496784, 0.0)
table["199Hg"] = (1, 16.87, 0.87621937, 4.8457916, 17.910822, 0.0)
table["203Ti"] = (1, 29.524, 2.80983305, 15.5393338, 57.123200, 0.0)
table["205Ti"] = (1, 70.476, 2.83747094, 15.6921808, 57.683838, 0.0)
#table["207Pb"] = (1, 22.1, 1.00906, 5.58046, 20.920599, 0.0)
#erroneous values in 2001 paper, corrected here
table["207Pb"] = (1, 22.1, 1.02638, 5.67625, 20.920599, 0.0)

# quadrupolar
table["2H"] = (2, 0.0115, 1.21260077, 4.10662791, 15.350609, 0.2860)
table["6Li"] = (2, 7.59, 1.1625637, 3.9371709, 14.716086, -0.0808)
table["7Li"] = (3, 92.41, 4.20407505, 10.3977013, 38.863797, -4.01)
table["9Be"] = (3, 100, -1.520136, -3.759666, 14.051813, 5.288)
table["10B"] = (6, 19.9, 2.0792055, 2.8746786, 10.743658, 8.459)
table["11B"] = (3, 80.1, 3.4710308, 8.5847044, 32.083974, 4.059)
table["14N"] = (2, 99.632, 0.57100428, 1.9337792, 7.226317, 2.044)
table["17O"] = (5, 0.038, -2.24077, -3.62808, 13.556457, -2.558)
table["21Ne"] = (3, 0.27, -0.854376, -2.11308, 7.894296, 10.155)
table["23Na"] = (3, 100, 2.8629811, 7.0808493, 26.451900, 10.4)
table["25Mg"] = (5, 10.00, -1.01220, -1.63887, 6.121635, 19.94)
table["27Al"] = (5, 100, 4.3086865, 6.9762715, 26.056859, 14.66)
table["33S"] = (3, 0.76, 0.8311696, 2.055685, 7.676000, -6.78)
table["35Cl"] = (3, 75.78, 1.061035, 2.624198, 9.797909, -8.165)
table["37Cl"] = (3, 24.22, 0.8831998, 2.184368, 8.155725, -6.435)
table["39K"] = (3, 93.2581, 0.50543376, 1.2500608, 4.666373, 5.85)
table["40K"] = (8, 0.0117, -1.4513203, 1.5542854, 5.802018, -7.3)
table["41K"] = (3, 6.7302, 0.27739609, 0.68606808, 2.561305, 7.11)
table["43Ca"] = (7, 0.135, -1.494067, -1.803069, 6.730029, -4.08)
table["45Sc"] = (7, 100, 5.3933489, 6.5087973, 24.291747, -22.0)
table["47Ti"] = (5, 7.44, -0.93294, -1.5105, 5.637534, 30.2)
table["49Ti"] = (7, 5.41, -1.25201, -1.51095, 5.639037, 24.7)
table["50V"] = (12, 0.250, 3.6137570, 2.6706490, 9.970309, 21.0)
table["51V"] = (7, 99.750, 5.8380835, 7.0455117, 26.302948, -5.2)
table["53Cr"] = (3, 9.501, -0.61263, -1.5152, 5.652496, -15.0)
table["55Mn"] = (5, 100, 4.1042437, 6.6452546, 24.789218, 33.0)
table["59Co"] = (7, 100, 5.247, 6.332, 23.727074, 42.0)
table["61Ni"] = (3, 1.1399, -0.96827, -2.3948, 8.936051, 16.2)
table["63Cu"] = (3, 69.17, 2.8754908, 7.1117890, 26.515473, -22.0)
table["65Cu"] = (3, 30.83, 3.07465, 7.60435, 28.403693, -20.4)
table["67Zn"] = (5, 4.10, 1.035556, 1.676688, 6.256803, 15.0)
table["69Ga"] = (3, 60.108, 2.603405, 6.438855, 24.001354, 17.1)
table["71Ga"] = (3, 39.892, 3.307871, 8.181171, 30.496704, 10.7)
table["73Ge"] = (9, 7.73, -0.9722881, -0.9360303, 3.488315, -19.6)
table["75As"] = (3, 100, 1.858354, 4.596163, 17.122614, 31.4)
table["79Br"] = (3, 50.69, 2.719351, 6.725616, 25.053980, 31.3)
table["81Br"] = (3, 49.31, 2.931283, 7.249776, 27.006518, 26.2)
table["83Kr"] = (9, 11.49, -1.07311, -1.03310, 3.847600, 25.9)
table["85Rb"] = (5, 72.17, 1.6013071, 2.5927050, 9.654943, 27.6)
table["87Rb"] = (3, 27.83, 3.552582, 8.786400, 32.720454, 13.35)
table["87Sr"] = (9, 7.00, -1.2090236, -1.1639376, 4.333822, 33.5)
table["91Zr"] = (5, 11.22, -1.54246, -2.49743, 9.296298, -17.6)
table["93Nb"] = (9, 100, 6.8217, 6.5674, 24.476170, -32.0)
table["95Mo"] = (5, 15.92, -1.082, -1.751, 6.516926, -2.2)
table["97Mo"] = (5, 9.55, -1.105, -1.788, 6.653695, 25.5)
table["99Tc"] = (9, 0, 6.281, 6.046, 22.508326, -12.9)
table["99Ru"] = (5, 12.76, -0.7588, -1.229, 4.605151, 7.9)
table["101Ru"] = (5, 17.06, -0.8505, -1.377, 5.161369, 45.7)
table["105Pd"] = (5, 22.33, -0.760, -1.23, 4.576100, 66.0)
table["113In"] = (9, 4.29, 6.1124, 5.8845, 21.865755, 79.9)
table["115In"] = (9, 95.71, 6.1256, 5.8972, 21.912629, 81.0)
table["121Sb"] = (5, 57.21, 3.9796, 6.4435, 23.930577, -36.0)
table["123Sb"] = (7, 42.79, 2.8912, 3.4892, 12.959217, -49.0)
table["127I"] = (5, 100, 3.328710, 5.389573, 20.007486, -71.0)
table["131Xe"] = (3, 21.18, 0.8931899, 2.209076, 8.243921, -11.4)
table["133Cs"] = (7, 100, 2.9277407, 3.5332539, 13.116142, -0.343)
table["135Ba"] = (3, 6.592, 1.08178, 2.67550, 9.934457, 16.0)
table["137Ba"] = (3, 11.232, 1.21013, 2.99295, 11.112928, 24.5)
table["138La"] = (10, 0.090, 4.068095, 3.557239, 13.194300, 45.0)
table["139La"] = (7, 99.910, 3.1556770, 3.8083318, 14.125641, 20.0)

table["177Hf"] = (7, 18.60, 0.8997, 1.086, 4.007, 336.5)
table["179Hf"] = (9, 13.62, -0.7085, -0.6821, 2.517, 379.3)
table["181Ta"] = (7, 99.988, 2.6879, 3.2438, 11.989600, 317.0)
table["185Re"] = (5, 37.40, 3.7710, 6.1057, 22.524600, 218.0)
table["187Re"] = (5, 62.60, 3.8096, 6.1682, 22.751600, 207.0)
table["189Os"] = (3, 16.15, 0.851970, 2.10713, 7.765400, 85.6)
table["191Ir"] = (3, 37.3, 0.1946, 0.4812, 1.718, 81.6)
table["193Ir"] = (3, 62.7, 0.2113, 0.5227, 1.871, 75.1)
table["197Au"] = (3, 100, 0.191271, 0.473060, 1.729, 54.7)
table["201Hg"] = (3, 13.18, -0.7232483, -1.788769, 6.611583, 38.6)
table["209Bi"] = (9, 100, 4.5444, 4.3750, 16.069288, -51.6)

# WARNING, a different setting for biomolecules was proposed in a IUPAC/IUPAB recommendation in 1998
#     Markley et al. Recommendations for the presentation of NMR structures of proteins and nucleic acids. Journal of Molecular Biology (1998) vol. 280 (5) pp. 933-952
# This previous recommendation was only mentionning 2D 13C, 15N and 31P, but was using  different references
# 
# Usage is to use 1998 recommendation for proteins and nucleic acids.
# these are enforced when mode = "UIPAB" or "biomolecule"
#
# the alternative references are cited in table 4 of the IUPAC-2001 recommendation, and uses DSS as 1H reference
# this makes +2.6645 ppm shift in 13C.
# this makes a +380.4434 ppm shift in 15N.

table4={}
table4["2H"] = (2, 0.0115, 1.21260077, 4.10662791, 15.350608, 0.2860)
table4["13C"] = (1, 1.07, 1.216613, 6.728284, 25.144953, 0.0)
table4["15N"] = (1, 0.368, -0.49049746, -2.71261804, 10.132912, 0.0)
table4["31P"] = (1, 100, 1.95999, 10.8394, 40.480864, 0.0)

def freq(spin="1H", H_freq=100.00, mode="standard"):
    """ returns the frequency of the given spin, in a field were the proton 1H is at frequencey H_freq values 
    if mode is "standard" or "IUPAC", (default) the standard values are used
    if mode is "biomolecule" or "IUPAB" values are taken from previous paper
        Markley et al. Recommendations for the presentation of NMR structures of proteins and nucleic acids. Journal of Molecular Biology (1998) vol. 280 (5) pp. 933-952
        This previous recommendation was only mentionning 2D 13C, 15N and 31P, but was using  different references
        They are given in table 4 of 2001 paper

    Usage is to use 1998 recommendation (mode = "biomolecule") for proteins and nucleic acids.
        this makes +2.6645 ppm shift in 13C.
        this makes a +380.4434 ppm shift in 15N.
    
    """
    if mode in ("biomolecule", "IUPAB"):
        try:
            (a,b,c,d,fq,f)=table4[spin]
        except KeyError:
            (a,b,c,d,fq,f)=table[spin]
    else:
        (a,b,c,d,fq,f)=table[spin]
    return abs(fq/100.0*H_freq)

def freqB(spin="1H", Bo=14.092):
    """ returns the frequency in MHz of the given spin, in a field of intensity Bo Tesla
    """
    from math import pi
    (a,b,c,gama,e,f)=table[spin]
    return abs(10*gama*Bo/(2*pi))

def receptivity(spinans, spinref="1H"):
    """returns the receptivity of a given spin, relative to spinref"""
    (spin,  naturalAbondance, magneticMoment, magnetogyricRatio, fq,  quadrupoleMoment) = table[spinref]
    receptref = magnetogyricRatio**3 * (spin*(spin+1)) * naturalAbondance
    (spin,  naturalAbondance, magneticMoment, magnetogyricRatio, fq,  quadrupoleMoment) = table[spinans]
    recept = magnetogyricRatio**3 * (spin*(spin+1)) * naturalAbondance
    return(abs(recept/receptref))

def report(H_freq=600.0):
    """ report the spin table at a given field
    """
    from math import pi
# compute Bo
    (spin,  naturalAbondance, magneticMoment, magnetogyricRatio, ff,  quadrupoleMoment) = table["1H"]
    Bo = H_freq *2* pi / (10* magnetogyricRatio)
# load freq table to sort and inverse it
    t = {}
    for i in table.keys():
        t[freq(i, H_freq)] = i
    l = t.keys()
    l.sort()
    l.reverse()
#    print l
    print("===   NPK Spin table  at Bo = %6.2f Tesla  ======================================"%(Bo))
    print("name   :  spin      Freq (\Xi) MHz    Q in fm2     Abundance  in %       receptivity")
    print("====================================================================================")
    for f in l:
        s = t[f]
        (spin,  naturalAbondance, magneticMoment, magnetogyricRatio, ff,  quadrupoleMoment) = table[s]
        f2 = freqB(s, Bo)
        if (spin%2 == 0):   # spin is integer
            print("%5s  :  spin %1i     Freq %8.4f    Q : %6.3f   Abundance : %7.3f%%  recept : %5.2e "% \
                (s, spin//2, freq(s,H_freq),quadrupoleMoment,naturalAbondance,receptivity(s)))
        else:
            if spin == 1:
                print("%5s  :  spin %1i/2   Freq %8.4f                 Abundance : %7.3f%%  recept : %5.2e"% \
                    (s,spin, freq(s,H_freq),naturalAbondance,receptivity(s)))
            else:
                print("%5s  :  spin %1i/2   Freq %8.4f    Q : %6.2f   Abundance : %7.3f%%  recept : %5.2e"% \
                    (s,spin, freq(s,H_freq),quadrupoleMoment,naturalAbondance,receptivity(s)))
#        if (abs(1- (f2/freq(s,H_freq))) > 0.01):
#            print f2, freq(s,H_freq), 100*abs(1- (f2/freq(s,H_freq)))
