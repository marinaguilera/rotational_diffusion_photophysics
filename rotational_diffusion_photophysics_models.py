import numpy as np
import rotational_diffusion_photophysics as rdp

'''
Usefull full model of experimentally realistic STARSS experiments are implemented
in this script.
'''

################################################################################
## Common

numerical_aperture = 1.4 # 1.4
refractive_index = 1.518 # 1.518
# numerical_aperture = 0.1 #
# refractive_index = 1 #
# numerical_aperture = 1.1 # reduced NA
# refractive_index = 1.518 # water at 405nm
lmax = 6

### Create the detectors
detXY = rdp.PolarizedDetection(polarization=['x', 'y'],
                               numerical_aperture=numerical_aperture,
                               refractive_index=refractive_index,
                               )

### rsEGFP2 Models
### Create the fluorophore, with the photophysics
rsEGFP2_4states = rdp.NegativeSwitcher(extinction_coeff_on=[5260, 51560],
                                        extinction_coeff_off=[22000, 60],
                                        wavelength=[405, 488],
                                        lifetime_on=1.6e-9,
                                        quantum_yield_on_fluo=0.35,
                                        quantum_yield_on_to_off=1.65e-2,
                                        )

rsEGFP2_8states = rdp.NegativeSwitcher(extinction_coeff_on=  [  5260, 51560],
                                        extinction_coeff_off=[ 22000,    60],
                                        wavelength=          [   405,   488],
                                        lifetime_on=1.6e-9,
                                        lifetime_off=20e-12,
                                        quantum_yield_on_to_off=1.65e-2,
                                        quantum_yield_off_to_on=0.33,
                                        quantum_yield_on_fluo=0.35,
                                        starting_populations=[1,0,0,0,0,0,0,0],
                                        deprotonation_time_off=5.1e-6,
                                        protonation_time_on=48e-6,
                                        nspecies=8,
                                        quantum_yield_trans_to_cis_anionic=0.0165,
                                        quantum_yield_cis_to_trans_neutral=0.33,
                                        )

rsEGFP2_9states = rdp.NegativeSwitcher(extinction_coeff_on=  [  5260, 51560],
                                        extinction_coeff_off=[ 22000,    60],
                                        wavelength=          [   405,   488],
                                        lifetime_on=1.6e-9,
                                        lifetime_off=20e-12,
                                        quantum_yield_on_to_off=1.65e-2,
                                        quantum_yield_off_to_on=0.33,
                                        quantum_yield_on_fluo=0.35,
                                        starting_populations=[1,0,0,0,0,0,0,0,0],
                                        deprotonation_time_off=5.1e-6,
                                        protonation_time_on=48e-6,
                                        nspecies=9,
                                        quantum_yield_trans_to_cis_anionic=0.0165,
                                        quantum_yield_cis_to_trans_neutral=0.33,
                                        )

# Diffusion model
def isotropic_diffusion(tau):
    return rdp.IsotropicDiffusion(diffusion_coefficient=1/(6*tau))


# Diffusion model
iso100us = rdp.IsotropicDiffusion(diffusion_coefficient=1/(6*100e-6))
iso1000us = rdp.IsotropicDiffusion(diffusion_coefficient=1/(6*1000e-6))
 
################################################################################
## STARSS1
'''
STARSS modality 1 and 2 can be simulated in a relatively easy way because only
one experiment is perfomed. So the pulse scheme is just one and the output of
the rdp class is already giving the experimental observable.
'''

### Modeled as starss1 experiments 2020/07/21
exc405X488C = rdp.ModulatedLasers(power_density=[176.7e3, 17.6e3/6],
                                  wavelength=[405, 488],
                                  polarization=['x', 'xy'],
                                  modulation=[[0,1,0],[1,1,1]],
                                  time_windows=[10e-3, 250e-9, 3e-3],
                                  time0=10e-3+250e-9,
                                  numerical_aperture=numerical_aperture,
                                  refractive_index=refractive_index,
                                  )

# STARSS modality 1
starss1 = rdp.System(illumination=exc405X488C,
                     fluorophore=rsEGFP2_8states,
                     diffusion=iso100us,
                     detection=detXY,
                     lmax=lmax)

################################################################################
## STARSS2
### Create the illumination scheme
exc488Xsimple = rdp.ModulatedLasers(power_density=[5000],
                              wavelength=[488],
                              polarization=['x'],
                              modulation=[[1]],
                              time_windows=[3e-3],
                              time0=0,
                              numerical_aperture=numerical_aperture,
                              refractive_index=refractive_index,
                              )

exc488X = rdp.ModulatedLasers(power_density=[114.8e3, 19.6e3/6],
                              wavelength=[405, 488],
                              polarization=['xy', 'x'],
                              modulation=[[0,1,0,0],[1,1,0,1]],
                              time_windows=[10e-3, 250e-9, 0.5e-3, 3e-3],
                              time0=10e-3+250e-9+0.5e-3,
                              numerical_aperture=numerical_aperture,
                              refractive_index=refractive_index,
                              )

### Full system
starss2 = rdp.System(illumination=exc488X,
                     fluorophore=rsEGFP2_4states,
                     diffusion=iso100us,
                     detection=detXY,
                     lmax=lmax)


################################################################################
## STARSS3
'''
STARSS modality 3 is a bit more complicated because every data point is computed
from a pulse scheme with a given time delay between two 405 pulses.
'''
### Pulse scheme with tunable delay between a pair of 405 pulses
def starss3_illumination(delay=5e-6, 
                         numerical_aperture=numerical_aperture,
                         refractive_index=refractive_index):
    illumination = rdp.ModulatedLasers(power_density=[7.12e6/6/2, 14e3/6],
                                       wavelength=[405, 488],
                                       polarization=['x', 'xy'],
                                       modulation=[[0,0,1,0,1,0,0],[1,0,0,0,0,0,1]],
                                       time_windows=[10e-3, 0.6e-3-delay, 50e-9, delay, 50e-9, .2e-3, 3e-3],
                                       time0=10e-3+.6e-3+50e-9+50e-9+.2e-3,
                                       numerical_aperture=numerical_aperture,
                                       refractive_index=refractive_index,
                                       )
    return illumination

def starss3_illumination_1pulse(delay=5e-6):
    illumination = starss3_illumination(delay)
    illumination.modulation = np.array([[0,0,0,0,1,0,0],[1,0,0,0,0,0,1]], dtype='float')
    return illumination

### Circularly polarized detector
detC = rdp.PolarizedDetection(polarization=['xy'],
                              numerical_aperture=numerical_aperture,
                              refractive_index=refractive_index,
                              )

### Full experiment
starss3 = rdp.System(illumination=starss3_illumination(50e-6),
                     fluorophore=rsEGFP2_8states,
                     diffusion=isotropic_diffusion(100e-6),
                     detection=detC,
                     )

### STARSS modality 3 signal as repeated pulse scheme with tunable delayed
def starss3_detector_signals(delays,
                             experiment=starss3,
                             ):

    s = np.zeros(np.size(delays))
    t = np.linspace(0,1e-3,1000)
    t = np.append(t, [10e-3])  # add a time with a long delay (background)
    integration_lims = [0, 1e-3]

    for i in np.arange(np.size(delays)):
        illumination = starss3_illumination(delays[i])
        illumination_1pulse = starss3_illumination_1pulse(delays[i])

        # signal with only one pulse
        experiment.illumination = illumination_1pulse
        signal_1pulse = experiment.detector_signals(t)

        # signal with two pulses
        experiment.illumination = illumination
        signal = experiment.detector_signals(t)

        signal_minus_bg_1pulse = signal_1pulse - signal_1pulse[0,-1]
        signal_minus_bg = signal - signal[0,-1]
        t_sel = np.logical_and(t>=integration_lims[0], t<=integration_lims[1])

        s[i] = np.sum(signal_minus_bg[0,t_sel]) / np.sum(signal_minus_bg_1pulse[0,t_sel])
    return s, experiment


if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    delay = [0.5e-6, 1e-6,10e-6,100e-6,200e-6,499e-6]
    s, exp = starss3_detector_signals(delay)
    # print(s)
    plt.figure(dpi=300)
    plt.plot(delay, s)
    plt.ylim((1.5,1.75))
    plt.xscale('log')
    plt.xlabel('Delay (s)')
    plt.ylabel('Normalized Counts')
    plt.show()

    # label = 'starss3_NA0p1_rsEGFP24states_550KWcm2_100us_'
    # plt.savefig(label+'normalized_counts')
    # np.savez_compressed(label+'data',s)

    # plt.figure()
    # plt.plot(exp._c[0].T)
    # plt.ylim((-1,1))

    # plt.figure()
    # plt.imshow(np.log10(np.abs(exp._M[1,:,:,0,0])))
    # plt.colorbar()

    # plt.figure()
    # plt.imshow(np.log10(np.abs(exp._M[1,1,0,:,:])))
    # plt.colorbar()

    # from plot_data_sphere import plot_data_sphere
    # c0_10 = exp._c[2,:,10]
    # grid = rdp.vec2grid(c0_10)
    # ax = plt.subplot(projection='3d')
    # plot_data_sphere(ax, grid.data)