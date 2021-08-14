#! /usr/bin/env python
import numpy as np
from numpy import *
#import matplotlib.pyplot as plt
import matplotlib
from matplotlib.pyplot import *
import xml.etree.cElementTree as ET
from scipy import interpolate
import scipy
from scipy.optimize import curve_fit
import copy
import json
import sys
import os

mode = 'error'
''' The valid argument would be python thermoCreator 'name of desired species' 'string of contributed groups'
    for example : 'C2H6' 'C-(C)(H)3:2'  '''
R = 8.314 #J/kmol.k
#R = 1.9872 #Cal/mol.k
T0 = 298.15
sigma = 1e-3
Tmid = 1000
def chemkin(p_cp_1 , p_H_1 , p_S_1, p_cp_0, p_H_0, p_S_0):
    line2 = ''
    end2 = '    2' 
    for i in range(0 , p_cp_1.size):
        if p_cp_1[i]>0:
            sep = ' '
            line2 += sep + str("{:.8E}".format(p_cp_1[i]))
        if p_cp_1[i]<0:
            sep = ''
            line2 += sep + str("{:.8E}".format(p_cp_1[i]))
    line2 += end2
    
    line3 = ''
    end3 = '    3'
    if p_H_1 >0:
	    sep = ' '
	    line3 += sep + str("{:.8E}".format(p_H_1))
    if p_H_1 <0:
	    sep= ''
	    line3 += sep + str("{:.8E}".format(p_H_1))
    if p_S_1 >0:
        sep = ' '
        line3 += sep + str("{:.8E}".format(p_S_1))
    if p_S_1 <0:
        sep= ''
        line3 += sep + str("{:.8E}".format(p_S_1))
    for i in range(0,3):
        if p_cp_0[i]>0:
            sep = ' '
            line3 += sep + str("{:.8E}".format(p_cp_0[i]))
        if p_cp_0[i]<0:
            sep = ''
            line3 += sep + str("{:.8E}".format(p_cp_0[i]))
    line3 += end3
    
    line4 = ''
    end4 = '                   4'
    for i in range(3,5):
        if p_cp_0[i]>0:
            sep = ' '
            line4 += sep + str("{:.8E}".format(p_cp_0[i]))
        if p_cp_0[i]<0:
            sep = ''
            line4 += sep + str("{:.8E}".format(p_cp_0[i]))
    if p_H_0 >0:
        sep = ' '
        line4 += sep + str("{:.8E}".format(p_H_0))
    if p_H_0 <0:
        sep= ''
        line4 += sep + str("{:.8E}".format(p_H_0))
    if p_S_0 >0:
        sep = ' '
        line4 += sep + str("{:.8E}".format(p_S_0))
    if p_S_0 <0:
        sep= ''
        line4 += sep + str("{:.8E}".format(p_S_0))
    line4 += end4
    print ('Polynomials in Chemkin format are as follow:')
    print (line2)
    print (line3)
    print (line4)
    return line2, line3, line4

def nasaStr( p_cp , p_H , p_S):
    erg = '['
    sep = ''
    for i in range(0 , p_cp.size):
      erg += sep + str("{:.8E}".format(p_cp[i]))
      sep = ', '
    erg += sep + str("{:.8E}".format(p_H))
    erg += sep + str("{:.8E}".format(p_S))
    erg += ']'
    return erg
  
def nasa(p_cp , p_H , p_S):   
    sep = ''
    for i in range(0 , p_cp.size):
      sep += str(p_cp[i])
      sep += ', '
    sep += str(p_H)
    sep += ', '
    sep += str(p_S)
    return sep

def trangeStr(tmin , tmax):
    erg = '['
    sep = ', '
    erg += str( tmin )
    erg += sep + str(tmax)
    erg += ']'
    return erg

def func_cp(T, a0, a1, a2, a3, a4):
  cp = R * (a0 + a1*T + a2*T*T + a3*T*T*T + a4*T*T*T*T)
  return cp 

def func_H(T, a0, a1, a2, a3, a4, a5):
  H = (a0 + (a1/2)*T + (a2/3)*T*T + (a3/4)*T*T*T + (a4/5)*T*T*T*T + (a5/T))*(R*T)
  return H

def func_S(T,a0, a1, a2, a3, a4, a6):
  S = R*(a0*log(T) + a1*T + (a2/2)*T*T + (a3/3)*T*T*T + (a4/4)*T*T*T*T + a6)
  return S       

def a6_from_S( S, a, T):
        a6 = (S/R) - (a[0]*log(T)) - (a[1]*T) - ((a[2]/2)*T*T) - ((a[3]/3)*T*T*T) - ((a[4]/4)*T*T*T*T)
        return a6
    
def a5_from_H( H, a, T ):
        a5 = (H/R) - (a[0]*T) - ((a[1]/2)*T*T) - ((a[2]/3)*T*T*T) - ((a[3]/4)*T*T*T*T) - ((a[4]/5)*T*T*T*T*T)
        return a5
    
def thermofit(int0_T,int1_T,int0_cp, int1_cp, Hf_298,  T_mid, S0 ):
    
    
    sigma0 =np.ones(len(int0_T))
    for i in [i for i,x in enumerate(int0_T) if x == T0]:
        sigma0[i] = sigma * 0.01
    sigma0[-1] = sigma #    * 0.1

    sigma1 = np.ones(len(int1_T))
    sigma1[0] = sigma #* 0.1

    guess = [1.0 , 1.0 , 1.0 , 0.0 , 0.0]  # initial guess for coeffs a0-a4
    
    params0_int, params0_covariance_int = scipy.optimize.curve_fit( func_cp,  int0_T , int0_cp, guess , maxfev=2000 , sigma=sigma0 )

    params1_int, params1_covariance_int = scipy.optimize.curve_fit( func_cp, int1_T , int1_cp, guess , maxfev=2000 , sigma=sigma1 )

    # ----------------------------
    # calculation of integration constant of enthalpy a5 for first and second interval.
    a5_0 = a5_from_H( Hf_298, params0_int, T0)

    H_mid = func_H(T_mid, params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4], a5_0)
    a5_1 = a5_from_H(H_mid, params1_int, T_mid)
    # ---------------------------
    #Finding a6 from the entropy formulation without having given values for S(T) using standard entropy at
    #298.15K for first interval and with entropy of ending point of last interval to calculate a6 for second interval  
    a6_0 = a6_from_S( S0, params0_int, T0)

    S_mid = func_S(T_mid, params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4], a6_0)
    a6_1 = a6_from_S(S_mid, params1_int, T_mid)
    

    
    chemkin(params1_int , a5_1, a6_1, params0_int , a5_0 , a6_0)
    p0str = nasaStr( params0_int , a5_0 , a6_0 )
    p1str = nasaStr( params1_int , a5_1, a6_1 )
    
    t0str = trangeStr( int0_T[0] , T_mid )
    t1str = trangeStr( T_mid , int1_T[-1] )

    n0str = '(NASA(' + t0str + ',\n'+ p0str + ')'
    n1str = 'NASA(' + t1str + ',\n'+ p1str + ')),'

    #print ("polynomials are:")
    #print (nameOfSpecies)
    #print (n0str)
    #print (n1str)
    #print (func_H(298.15 , params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4], a5_0 ))
    a = func_cp(Tmid , params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4] ) - func_cp(Tmid , params1_int[0] , params1_int[1] , params1_int[2], params1_int[3] , params1_int[4] )
    #------------------
    b = func_H(Tmid, params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4],a5_0 ) - func_H(Tmid, params1_int[0] , params1_int[1] , params1_int[2], params1_int[3] , params1_int[4], a5_1 )
    #------------------
    c = func_S(Tmid, params0_int[0] , params0_int[1] , params0_int[2], params0_int[3] , params0_int[4], a6_0) - func_S(Tmid, params1_int[0] , params1_int[1] , params1_int[2], params1_int[3] , params1_int[4], a6_1)
    print ("Error at Tmid : ")
    print ("Cp / H  / S = " + \
      str(a) + " / " + \
      str(b) + " / "   + \
      str(c))
    
    nasa0 = nasa( params0_int , a5_0 , a6_0 )
    nasa1 = nasa( params1_int , a5_1, a6_1 )
    return nasa0, nasa1
    
    
def main_run():
    #print ('start of the code')
    nameOfSpecies = sys.argv[1]
    groupString = sys.argv[2]
    setupFile = 'groups_data.xml'
    eleRoot = ET.parse(setupFile)
    speciesData = eleRoot.find('speciesData')
    species = speciesData.findall('species')
    groupdatabase = []
    for s in species:
        s_name=s.attrib['name']
        groupdatabase.append(s_name)
    
    print ('string given: ' + str(groupString))
    print ('for this species: ' + str(nameOfSpecies))
    gls = groupString.split(' ')
    gl = [x.strip() for x in gls]
    listOfGroup = {}

    for thisgl in gl :
        smult = thisgl.split(':')
        smult = [x.strip() for x in smult]
        if (len(smult) == 2):
            groupName = smult[0]
            if (groupName in groupdatabase):
                groupValue = int(smult[1])
                listOfGroup[groupName] = groupValue
            else:
                raise Exception("Such a group does not exist in the database: " + groupName)  
        else :
            raise Exception("Uninterpreted Element")
            #print ('Uninterpreted Element: ' + str(thisgl))
    species_list = list(listOfGroup.keys())
    cp_0 = {}
    cp_1 = {}
    hf = {}
    s0 = {}
    T_0 = np.arange(200,1000+1, 100)
    T_1 = np.arange(1000,4000+1, 100)
    groupArray = ''
    #print (species)
    for i in range(len(species_list)):
        groupArray += species_list[i]
        groupArray += ':'
        groupArray += str(listOfGroup[species_list[i]])
        groupArray += ' '
    
    for s in species:
        s_name = s.attrib['name']   
        if (s_name in species_list ):
            thermo = s.find('thermo')
            nasa = thermo.findall('NASA')
            cp_0[s_name] = []   
            cp_1[s_name] = []
            
            for n in nasa:
                if (float(n.attrib['Tmax']) == Tmid):
                    poly = n.find('floatArray')
                    fitted_int0 = poly.text.split(',')
                    fitted_int0 = [float(i) for i in fitted_int0] 
                    for temp in T_0 :
                        cp_0[s_name].append(listOfGroup[s_name]* func_cp(temp,fitted_int0[0],fitted_int0[1],fitted_int0[2],fitted_int0[3],fitted_int0[4]))
                    hf_298 = func_H(T0,fitted_int0[0],fitted_int0[1],fitted_int0[2],fitted_int0[3],fitted_int0[4],fitted_int0[5] ) 
                    hf[s_name] = listOfGroup[s_name] * hf_298 
                    s0_298 = func_S(T0,fitted_int0[0],fitted_int0[1],fitted_int0[2],fitted_int0[3],fitted_int0[4],fitted_int0[6] ) 
                    s0[s_name]  = listOfGroup[s_name] * s0_298
                elif (float(n.attrib['Tmin']) == Tmid):
                    poly = n.find('floatArray')
                    fitted_int1 = poly.text.split(',')
                    fitted_int1 = [float(i) for i in fitted_int1]
                    for temp in T_1 :
                        cp_1[s_name].append(listOfGroup[s_name] * func_cp(temp,fitted_int1[0],fitted_int1[1],fitted_int1[2],fitted_int1[3],fitted_int1[4]))   
        
    erg_hf = sum(list(hf.values()))
    erg_s = sum(list(s0.values()))
    array_0 = np.asarray(list(cp_0.values()))
    erg_cp_0 = array_0.sum(axis=0)
    array_1 = np.asarray(list(cp_1.values()))
    erg_cp_1 = array_1.sum(axis=0)
    
    s_data = {}
    s_data['name'] = nameOfSpecies
    s_data['cp_0'] = list(erg_cp_0)
    s_data['cp_1'] = list(erg_cp_1)    
    s_data['Hf_298'] = erg_hf
    s_data['S0'] = erg_s

    print ('Standard enthalpy of formation: ' + str(s_data['Hf_298']/1000) )
    print ('Standard entropy: ' + str(s_data['S0']))
    nasa_poly = thermofit(T_0,T_1,s_data['cp_0'] , s_data['cp_1'] ,s_data['Hf_298'] ,  Tmid, s_data['S0'] )
    
''' call the main function only within the main Thread'''
if __name__ == "__main__":
   main_run()
