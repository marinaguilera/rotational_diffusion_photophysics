import rotational_diffusion_photophysics as rdp
from os.path import dirname, join as pjoin
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

mat_fileName =(r'C:\Users\Guillem\Documents\Work\Triplet Studies\Results\Beads library\mat_files\2023_07_03\eosinY_N2_p008_30nm_beads_processed_bkg_free.mat')
mat_contents = sio.loadmat(mat_fileName)

global t
global r
global i_par
global i_perp
global system

t = mat_contents['t'][0:1000,0] * 1e-6 # [s]
r = mat_contents['r'][0:1000,0]
i_par = mat_contents['i_par'][0:1000,0]
i_perp = mat_contents ['i_perp'][0:1000,0]
# tr = np.zeros((4, t.size))
norm = np.max(i_par)
i_par = i_par/norm
i_perp = i_perp/norm

def log_bin_data(data):
    data_bin = np.zeros((data.size))
    power = 0
    position = 0
    binsize = 1

    while position+binsize+1 < data.size:
        data_bin[power] = np.mean(data[int(position):(int(position)+int(binsize))])
        power = power+1
        position = position+binsize
        binsize = np.round(np.float_power(1.1,power))

    data_bin = data_bin[0:power-1]
    return data_bin

# Intialize empty experiment
lmax = 6
system = rdp.System(illumination=[],
                    fluorophore=[],
                    diffusion=[],
                    detection=[],
                    lmax= lmax)

def fitting_model(t,k):
    # Define model ingridients
    ingridients = model_ingridients(k)
    system.illumination = ingridients[0]
    system.detection = ingridients[1]
    system.fluorophore = ingridients[2]
    system.diffusion = ingridients[3]

    # Calculate signals
    s = system.detector_signals(t)
    s = s -s[:,-1][:,None] # remove background
    r_theoretical = rdp.anisotropy(s)
    i_par_theo = s[0]
    i_perp_theo = s[1]

    # Find the best variable projection amplitude for the intensities in each channel
    A = var_proj_amplitude(i_par_theo, i_perp_theo, i_par, i_perp)
    i_par_theo = i_par_theo / A[0]
    i_perp_theo = i_perp_theo / A[0]

    model_output = np.concatenate((r_theoretical, i_par_theo, i_perp_theo))
    return model_output

def var_proj_amplitude(i_par_theo, i_perp_theo, i_par, i_perp):
    E = np.zeros((2,np.size(i_par)))
    T = np.zeros((2,np.size(i_par)))

    E[0] = i_par
    E[1] = i_perp
    T[0] = i_par_theo
    T[1] = i_perp_theo

    E = np.array([E.flatten()])
    T = np.array([T.flatten()])

    E_inv = np.linalg.pinv(E)
    A = T.dot(E_inv)
    return A

def model_ingridients(k):
    # Objective parameters
    numerical_aperture = 1.4
    refractive_index = 1.518

    # Illumination scheme
    exc488X = rdp.ModulatedLasers(power_density=[16e3*0.66],
                                  wavelength=[488],
                                  polarization=['x'],
                                  modulation=[[1]],
                                  time_windows=[1e-3],
                                  time0=0,
                                  numerical_aperture=numerical_aperture,
                                  refractive_index=refractive_index)
    # Detectors
    detXY = rdp.PolarizedDetection(polarization=['x','y'],
                                   numerical_aperture=numerical_aperture,
                                   refractive_index=refractive_index)

    # Fluorophore
    eosinY_4states = rdp.TripletProbe(extincion_coeff_exc=[30000],
                                      wavelength=[488],
                                      lifetime_exc1=2.1e-9,
                                      lifetime_dark1=40e-6,
                                      lifetime_dark2=50,
                                      quantum_yield_fluo=0.57,
                                      quantum_yield_isc= 0.07,
                                      quantum_yield_bleach=2e-2,
                                      starting_populations=[1,0,0,0])

    # Diffusion model
    isodiff = rdp.IsotropicDiffusion(diffusion_coefficient=1/(6*k))
    return exc488X, detXY, eosinY_4states, isodiff

def model_ingridients_triplet_probe(k):

    # Objective parameters
    numerical_aperture = 1.4
    refractive_index = 1.518

    # Illumination scheme
    exc488X = rdp.ModulatedLasers(power_density=[16e3*0.66],
                                  wavelength=[488],
                                  polarization=['x'],
                                  modulation=[[1]],
                                  time_windows=[1e-3],
                                  time0=0,
                                  numerical_aperture=numerical_aperture,
                                  refractive_index=refractive_index)

    # Detectors
    detXY = rdp.PolarizedDetection(polarization=['x','y'],
                                   numerical_aperture=numerical_aperture,
                                   refractive_index=refractive_index)

    # Fluorophore
    eosinY_4states = rdp.TripletProbe(extincion_coeff_exc=[30000],
                                      wavelength=[488],
                                      lifetime_exc1=2.1e-9,
                                      lifetime_dark1=28e-6,
                                      lifetime_dark2=50,
                                      quantum_yield_fluo= .57,
                                      quantum_yield_isc= 0.085,
                                      quantum_yield_bleach= k[0],
                                      starting_populations=[1,0,0,0])
    # Diffusion model
    isodiff = rdp.IsotropicDiffusion(diffusion_coefficient=1/(6*16e-3))
    return exc488X, detXY, eosinY_4states, isodiff

def function_2_fit(t,k):
    r_m, i_par_m, i_perp_m = fitting_model(t,k)

    # Log-bin data
    r_bin = log_bin_data(r_m)
    i_par_bin = log_bin_data(i_par_m)
    i_perp_bin = log_bin_data(i_perp_m)

    model_output = np.concatenate((r_bin, i_par_bin, i_perp_bin))
    return model_output

bias = 1
t_max = 100e-6
k0 = [3.5e-6]

t_sel = t< t_max
t_sim = t[t_sel]

r = r[t_sel]
i_par = i_par[t_sel]
i_perp = i_perp[t_sel]

r = log_bin_data(r)
i_par = log_bin_data(i_par)
i_perp = log_bin_data(i_perp)

t_fit = log_bin_data(t_sim)

ydata = np.concatenate((r, i_par, i_perp))

fit_output = curve_fit(fitting_model, xdata=t_fit, ydata=ydata, p0=k0, method='trf')

k = fit_output['popt']

print(k)










