import sys
sys.path.insert(0, './')
# sys.path.insert(0, '../')
import rotational_diffusion_photophysics as rdp
import numpy as np
import matplotlib.pyplot as plt
import codetiming

# Maximum l value for SH expansion
lmax = 6  # 6 is enough for good accuracy, with reasonable level of saturation

# Objective parameters
numerical_aperture = 1.4
refractive_index = 1.518

# Create the illumination scheme
exc488X = rdp.SingleLaser(power_density=800,
                         wavelength=488,
                         polarization='x',
                         numerical_aperture=numerical_aperture,
                         refractive_index=refractive_index,
                         )

# Create the detectors
detXY = rdp.PolarizedDetection(polarization=['x', 'y'],
                               numerical_aperture=numerical_aperture,
                               refractive_index=refractive_index,
                               )

# Create the fluorophore, with the photophysics
rsEGFP2_100us = rdp.NegativeSwitcher(cross_section_on_blue=2e-18,
                                     lifetime_on=3e-9,
                                     quantum_yield_on_to_off=0.001,
                                     diffusion_coefficient=1/(6*200e-6))

# Full system and experiment
starss2 = rdp.System(illumination=exc488X,
                     fluorophore=rsEGFP2_100us,
                     detection=detXY,
                     lmax=lmax)

# Visualize the photon fluexes wigner3j prod SH coefficients
F = starss2._F
plt.figure()
plt.imshow(F[0])

# Compute the diffusion and kinetics matrix
M = starss2.diffusion_kinetics_matrix()
plt.figure()
plt.imshow(M[0,0])

# Show the SH coefficients for the collection functions of the detectors
c_det = starss2._c_det
plt.figure()
plt.plot(c_det.T)

# Solve the time evolution and time it for pesrformances
t = np.logspace(-6, -3.5, 32)
tic = codetiming.Timer()
tic.start()
c = starss2.solve(t)
tic.stop()

# Compute signals, it includes the analytical solution of the problem
tic.start()
s = starss2.detector_signals(t)
tic.stop()

plt.figure()
plt.plot(t, c[1].T)
plt.xlabel('Time (s)')
plt.xscale('log')

plt.figure()
plt.plot(t, s.T)
plt.ylabel('Fluorescence Signal (a.u.)')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.legend(['Parallel', 'Perpendicular'])

plt.figure()
plt.plot(t, rdp.anisotropy(s))
plt.xlabel('Time (s)')
plt.xscale('log')
plt.ylabel('Anisotropy')
plt.ylim([0,0.4])

plt.figure()
plt.plot(t, rdp.anisotropy(s))
plt.xlabel('Time (s)')
plt.yscale('log')
plt.ylabel('Anisotropy')

print(rdp.anisotropy(s)[-1])

plt.show()
