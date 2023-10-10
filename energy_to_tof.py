# -*- coding: utf-8 -*-
"""
Convert Neutron Energy Spectrum to Time-of-Flight Spectrum

Created on Mon Oct  9 14:54:56 2023

@author: ziyuz
"""

import numpy as np
import matplotlib.pyplot as plt
import dress
mn = 1.674e-27 # mass of neutron in kg

#%%
reaction = dress.reactions.DTNHe4Reaction()
spec_calc = dress.SpectrumCalculator(reaction, n_samples=1e7)

# Distribution of 1st reactant (aka reactant a)
# -----------------------------------------
particle_a = spec_calc.reactant_a.particle
Ta = 3.0              # temperature of reactant a (keV)
dist_a = dress.dists.MaxwellianDistribution(Ta, particle_a)

# Distribution of 2nd reactant (aka reactant b)
# -----------------------------------------
particle_b = spec_calc.reactant_b.particle
Tb = 3.0
dist_b = dress.dists.MaxwellianDistribution(Tb, particle_b)

# Sample velocities from the distributions
# ------------------ ----------------------

spec_calc.reactant_a.v = dist_a.sample(spec_calc.n_samples)
spec_calc.reactant_b.v = dist_b.sample(spec_calc.n_samples)

# Calculate spectrum in a given direction
# ---------------------------------------
spec_calc.u = [0,0,1]       # Tip: try playing with the emission direction to see how this impacts the spectrum
spectrum, bin_centers = spec_calc()

plt.plot(bin_centers, spectrum)
plt.xlabel('Neutron energy (keV)')
plt.ylabel('Spectrum (au)')

#%%
#ntof
dist = 20.0 #vary distance here
nToF = np.zeros(len(bin_centers))
for i in range(len(bin_centers)):
    nToF[i] = dist/np.sqrt(2*bin_centers[i]*1.6022e-16/mn)*1e9 #in nanoseconds
plt.plot(nToF, spectrum)
plt.xlabel('Time of Flight (ns)')
plt.ylabel('Spectrum (au)')
plt.legend([dist])

#%%
#estimate FWHM assuming Gaussian
mean = np.average(bin_centers,weights=spectrum)
std = np.sqrt(np.average((bin_centers-mean)**2,weights=spectrum))
print(2.355*std)