#First attempt at a header programme for the set-up and comparison of the two
#models



#Initialisation
import numpy as np
import scipy
import math
import matplotlib
import matplotlib.pyplot as plt
#To track figure files
figno = 0
#Set up the required vector arrays
variables = []
measures = []
immutables = []
parameters = []

#first set the constants etc.
#ie
#set up a vector of the relevant values, print and save this
#(Assuming we can the import the values back into the fortran programms
#which will require a rewrite)

#PARAMETERS

#SET PARAMETERS THAT YOU WILL VARY HERE
#Such as

#The Matter Omega Parameter
omega_m = 0.3
#The Dark Energy Omega Parameter
omega_l = 1 - omega_m
#The Present Day Hubble Constant
h0 = 69

variables = [omega_m,omega_l,h0]

print('variables', variables)

#AND JUST COMPARE PHYSICAL CONSTANTS FOR NOW

#NUMERICAL PARAMETERS
#TIME steps

#INITIAL Conditions

#NO START WITH UNITS...

#unit measure of space
lu = 1
#Time unit
tu = 1
#mass unit
mu = 1

measures = [lu,tu,mu]

print('measures', measures)

pi = np.pi

print('the elements of the array should be (pi,)')

print('pi, defined by numpy.pi is ', pi)

lightspeed = 2.98*(mu/tu)

immutables = [lightspeed]

print('immutables', immutables)



parameters_t1 = [variables,measures,immutables]
print('parameters_t1:',parameters_t1)

parameternames_t1 = ['variables','measures','imutables']

prmtrnms_t2 = np.array(['omega_m','omega_l','h0','lu','tu','mu','ls(c)'])

print('names (t2)',prmtrnms_t2)

parameters_t2 = np.concatenate([variables,measures,immutables])
print('parameters_t2 (concatenated):',parameters_t2)

parameters_t3 = np.array(variables+measures+immutables)
print('parameters_t3 (added):',parameters_t3)

szpn = prmtrnms_t2.shape
szprmtr = parameters_t3.shape

print('shape of prmtr name array:', szpn)
print('shape of prtmr value array:', szprmtr)

prmtrdata = [prmtrnms_t2, parameters_t3]

#np.savetxt('parameters.txt',prmtrdata,delimiter=',',fmt='%d')
#np.savetxt('parameters.txt',prmtrdata,delimiter=',')
np.savetxt('parameterst123.txt',parameters_t3,delimiter=',',fmt='%d')

prmtrdt_rslt = str(prmtrdata)

f = open('parameters_t1.txt','w+')
f.write(prmtrdt_rslt)
f = close()

#f = open('parameters_t2.txt','w+')
#f.write(prmtrdata)
#f = close()

#np.save('paramter_data.txt', prmtrdata)


#COMPILE AND RUN FORTRAN PROGRAMMES

#EXTRACT RELEVANT SAVED DATA

#COMPARE RESULTS

#PRINT AND SAVE

#END
