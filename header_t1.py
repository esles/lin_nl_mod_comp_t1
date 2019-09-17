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


#SET PARAMETERS THAT YOU WILL VARY HERE
#SUCH AS

#The Matter Omega Parameter
omega_m = 0.3
#The Dark Energy Omega Parameter
omega_l = 1 - omega_m
#The Present Day Hubble Constant
h0 = 69

variables = [omega_m,omega_l,h0]

#AND JUST COMPARE PHYSICAL CONSTANTS FOR NOW


#NO START WITH UNITS...

#unit measure of space
lu = 1
#Time unit
tu = 1
#mass unit
mu = 1

measures = [lu,tu,mu]

pi = np.pi

print('the elements of the array should be (pi,)')

print('pi, defined by numpy.pi is ', pi)

lightspeed = 2.98*(mu/tu)

immutables = [lightspeed]

parameters = [variables,measures,immutables]
