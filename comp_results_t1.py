#
#(Python) Programme to process (ie input) and present (ie plot and compare) the
#data produced by the (Fortran) programmes lin_mod_v1.f90 and simsilsun_v1.f90
#written by me and Krzysztof Bolejko respectively (with modifications by me)
#

#
#The goal is to eventually merge this into a (python) header programme that can
#set the parameters, compile the fortran programmes, run the script and plot
#the results in one go.
#But we'll worry about this later


#Initialisation

#Import Libraries
import numpy as np
import matplotlib.pyplot as plt


#To track figure files
figno = 1


#The names are only temporary.
#The names are hardwired (ie need to know in advance what files exist in the
#current directory to import)

simsilun_density = np.loadtxt(fname='density')
simsilun_initial_density = np.loadtxt(fname='initial_density')

#Debug/Tackers
dtyp = simsilun_density.dtype
dsh = simsilun_density.shape

#Debug/Tackers
print('Misc. Info.')
print(dtyp)
print(dsh)

density_c0 = simsilun_density[:,0]
density_c1 = simsilun_density[:,1]
density_c2 = simsilun_density[:,2]

initial_density_c0 = simsilun_initial_density[:,0]
initial_density_c1 = simsilun_initial_density[:,1]

print('producing figure',figno)
#Initialise the figure (frame?)
plt.figure(figno)

#Add the data
plt.plot(density_c0,initial_density_c1, label='Initial')
plt.plot(density_c0,density_c1, label = 'Silent')
plt.plot(density_c0,density_c2, label = 'EdS')

#Add Title
plt.title('title')
#Add Axis Labels
plt.xlabel('x-axis: R')
plt.ylabel('y-axis: $\delta$')
#Add Legend
plt.legend(loc = 'lower right')


#Show the image
print('showing figure',figno)
plt.show()

#Save the Figure
#plt.savefig('test_comp.eps')


#Increment Figure Counter
figno = figno + 1
