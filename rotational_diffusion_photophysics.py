import numpy as np
import pyshtools as sht
import spherical as sf  # currently not used, might be removed in the future
import matplotlib.pyplot as plt  # used only for plotting tools

'''
This module solves analytically the rotational diffusion and kinetics
in flurescence experiments with complex photophysics.
An arbitrary kinetic scheme can be implemented using fluorophore classes.
The illumination of the experiment is implemented as a pulse scheme with
square wave modulations.
Detection and illumination can be done in arbitrary NA objectives, using 
the Axelrod correction.
'''

################################################################################
# Classes for programming experiments
################################################################################

class System:
    def __init__(self,
                 fluorophore,
                 diffusion,
                 illumination,
                 detection=[],
                 lmax=6,
                 norm='4pi',
                 ):
        # The simulation requires the expansion of all angular functions in 
        # spherical harmonics (SH). lmax is the cutoff for the maximum l quantum
        # number of the SH basis set, i.e. 0 <= l <= lmax.
        #
        # The '4pi' normalization is reccomended because the coefficients for
        # l=0 and m=0 for the angular probability distributions give directly
        # the state populations.
        self.lmax = lmax
        self.norm = norm

        # Compute arrays with quantum numbers l and m of the SH basis set.
        self._l, self._m = quantum_numbers(self.lmax)

        # Compute wigner3j coefficients. These are necessary for the evaluation
        # of the product of angular functions from the SH expansion 
        # coefficients.
        self._wigner3j_prod_coeffs = wigner_3j_prod_3darray(self._l, self._m)

        # Import the classes containing the parametrization and characteristics
        # of fluorophore, diffusion model, illumination, and detection.
        self.fluorophore = fluorophore
        self.diffusion = diffusion
        self.illumination = illumination
        self.detection = detection

        '''
        Small description of hidden variables. 
        These are used mainly for debugging.
        self._F is the photon flux prod matrix, when multiplied by the cross
            section gives the absorption rate prod matrix as a function of the 
            angle.
        self._c_exc are SH coefficients for the orientational functions behind F. 
            It can be interpred as the number of photons that would be absorbed 
            by a molecule with cross section 1cm2 as a function of the 
            orientation of the dipole.
        self._D Diffusion matrix rate model for each species.
        self._K Kinetic matrix rate model for each time window and each pair 
            of species.
        self._M = M Full diffusion-kinetics matrix for each time window and each 
            pair of species.
        self._c0 SH coeffs of starting orientational populations
        self._c SH coeffs of populations computed at the derired times
        self._c_det SH coeffs of the detectors collection functions
        self._s fluorescent signal for each detector
        '''
        return None

    def detector_signals(self, time):
        # Compute the flurescence signal registered in the detectors
        # Firt, solve the diffusion/kinetics problem and get the SH coefficients
        # of populations.
        c = self.solve(time)

        # Multiply the populations by the quantum yields of fluorescence.
        # Usually, only the fluorescent state will have a quantum yield 
        # different from zero.
        c_fluo = c * self.fluorophore.quantum_yield_fluo[:,None,None]

        # Compute the SH coefficients for the collection functions of 
        # the detectors
        c_det = self.detection.detector_coeffs(self._l, self._m)
        ndetectors = c_det.shape[0]

        # Initialize signal array and compute signals
        s = np.zeros( (ndetectors, time.size) )
        for i in range(ndetectors):
            s[i] = np.sum(c_det[i][None,:,None] * c_fluo, axis=(0,1))
        if self.norm == '4pi':
            s = s/4*np.pi

        # Save variables, mainly for debugging.
        self._c_det = c_det
        self._s = s
        return s
            
    def solve(self, time):
        # Solve the time evolution of the system including the pulse scheme.
        
        # Compute the diffusion kinetics matrix
        # Matrix M contains all the information about the time evolution of the 
        # system in every time window.
        M = self.diffusion_kinetics_matrix()

        # Get the initial conditions for the experiment from the fluorophore
        # class.
        c0 = self.fluorophore.starting_coeffs(self._l, self._m)

        # Optimization, remove all odd l value coefficients from c0 and M.
        # We will need to add them back, after the solution.
        # This semplification is usefull when only light matter interaction can
        # affect the orientational probablities.None
        # print(M.shape)
        # M, c0 = remove_odd_coeffs(M, c0, self._l)
        # print(M.shape)

        # Shift the time with time0, so it is represented in the laboratory
        # time (starting with the time zero of pulse sequence).
        time_lab = time + self.illumination.time0
        time_mod = np.cumsum(self.illumination.time_windows)

        # Insert a zero at the beginning of the modulation time sequence and 
        # extend last time window to fit the whole requested times to sovle.
        time_mod = np.insert(time_mod, 0, 0 )
        time_mod[-1] = np.max(time_lab)

        # Solve evolution in every time window
        c = np.zeros((self.fluorophore.nspecies,
                      self._l.size,
                      time.size))
        for i in np.arange(self.illumination.nwindows):
            # Selection of time axis inside the current time window
            time_sel = np.logical_and(time_lab >= time_mod[i],
                                      time_lab <= time_mod[i+1])

            # Select the time points inside the current time window
            time_i = time_lab[time_sel]

            # Add to the time points the end of the time window
            # so we can comute the starting coefficients for the next window.
            time_i = np.append(time_i, time_mod[i+1])

            # Shift the local time axis so the zero concides with the beginning
            # of the current time window
            time_i = time_i - time_mod[i]

            # Solve the time evolution
            c_i, _, _ = solve_evolution(M[i], c0, time_i)

            # Save results and update initial conditions for the next window
            c[:,:,time_sel] = c_i[:,:,:-1]
            c0 = c_i[:,:,-1]
        
        # Add back zeros in place of odd l coeffs
        # M, c0, c = add_odd_coeffs_zeros(M, c0, c, self._l)

        # Save variables, mainly for debugging.
        self._c0 = c0
        self._c = c
        return c

    def diffusion_kinetics_matrix(self):
        # Preliminary computations based on the illumination class.
        # Photon flux product coefficients based on wigner3j symbols
        F, c_exc = self.illumination.photon_flux_prod_coeffs(
                                                    self._l,
                                                    self._m,
                                                    self._wigner3j_prod_coeffs)

        # Compute the rotational diffusion matrix, it will not change with the
        # time windows of laser modulation.
        # The diffusion model is provided by the fluorophore class.
        # In the future, if more complex diffusion models will be implemented,
        # it might be convinient to create a separate diffusion class.
        D = self.diffusion.diffusion_matrix(self._l, self._m,
                self.fluorophore.nspecies)

        # Prepare the rotational diffusion matrix for each time window
        nwindows = self.illumination.nwindows
        nspecies = self.fluorophore.nspecies
        M = np.zeros( (nwindows,
                       nspecies, nspecies,
                       self._l.size, self._l.size) )
        K = np.zeros( (nwindows,
                       nspecies, nspecies,
                       self._l.size, self._l.size) )
        for i in np.arange(nwindows):
            K[i] = self.fluorophore.kinetics_matrix(self._l, self._m,
                        F * self.illumination.modulation[:,i][:,None,None],
                        self.illumination.wavelength)
            M[i] = diffusion_kinetics_matrix(D, K[i])

        # Save variables, mainly for debugging.
        self._c_exc = c_exc
        self._F = F
        self._D = D
        self._K = K
        self._M = M
        return M

def anisotropy(signals):
    return (signals[0] - signals[1]) / (signals[0] + 2*signals[1])

##########################
# Diffusion Models Classes
##########################
class IsotropicDiffusion:
    def __init__(self,
                 diffusion_coefficient=14e-9,  # GFP rotational diffusion
                 ):
        # Rotational diffusion coefficient in Hertz
        self.diffusion_coefficient = diffusion_coefficient  # [Hz]

    def diffusion_matrix(self, l, m, nspecies):
        # Compute the diffusion matrix
        # Here an isotropic diffusion model is employed, every state has also
        # the same rotational diffusion properties.
        D = isotropic_diffusion_matrix(l, m,
                self.diffusion_coefficient,
                nspecies)
        return D

#####################
# Fluorophore classes
#####################
class NegativeSwitcher:
    def __init__(self,
                 extinction_coeff_on=[0, 0],
                 extinction_coeff_off=[0, 0],
                 wavelength=[405, 488],
                 lifetime_on=3e-9,
                 lifetime_off=16e-12,
                 quantum_yield_on_to_off=0.001,
                 quantum_yield_off_to_on=0.2,
                 quantum_yield_on_fluo=1,
                 starting_populations=[1,0,0,0],
                 deprotonation_time_off = 15e-6,  # for 6 states model
                 protonation_time_on = 150e-6,  # for 6 states model
                 quantum_yield_trans_to_cis_anionic=0,  # for 8 states model
                 quantum_yield_cis_to_trans_neutral=0, # for 8 states model
                 nspecies=4):
        # Cross section in cm2 of absorptions
        epsilon2sigma = 3.825e-21  # [Tkachenko2007, page 5]
        self.extinction_coeff_on = np.array(extinction_coeff_on)  # [M-1 cm-1]
        self.extinction_coeff_off = np.array(extinction_coeff_off)  # [M-1 cm-1]
        self.cross_section_on = self.extinction_coeff_on * epsilon2sigma  # [cm2]
        self.cross_section_off = self.extinction_coeff_off * epsilon2sigma  # [cm2]
        self.wavelength = np.array(wavelength)  # [nm]

        # Lifetime of the on excited state in seconds
        # Assumption: lifetime on and off are the same for the same protonation
        # state. This might not be the case, especially because cis_neutral
        # species is not fluorescent, so most likely it will have a shorter
        # lifetime. For small enough excitation kinetic rates, if we don't have
        # accumulation of exited state species, then it doesn't matter much.
        self.lifetime_on = lifetime_on  # [s]
        self.lifetime_off = lifetime_off  # [s]

        # Quantum yield of an off-switching event from the on excited state and
        # fluorescence from the on state
        self.quantum_yield_on_fluo = quantum_yield_on_fluo
        self.quantum_yield_on_to_off = quantum_yield_on_to_off  # cis_to_trans_anionic 
        self.quantum_yield_off_to_on = quantum_yield_off_to_on  # trans_to_cis_neutra

        # Quantum yeilds and Protonation and deprotonation times for 6 and 8
        # states models.
        self.quantum_yield_cis_to_trans_neutral = quantum_yield_cis_to_trans_neutral
        self.quantum_yield_trans_to_cis_anionic = quantum_yield_trans_to_cis_anionic
        self.protonation_time_on = protonation_time_on
        self.deprotonation_time_off = deprotonation_time_off

        # Label describing the fluorophore type
        # Number of states in the kinetic model
        self.nspecies = nspecies

        # Index of the fluorescent state
        # Here, only one fluorescent state is assumed.
        # This is not necessary in general and could be extended in the future.
        self.quantum_yield_fluo = np.zeros((nspecies))
        self.quantum_yield_fluo[1] = self.quantum_yield_on_fluo

        # Population at the beginning of the experiment
        self.starting_populations = starting_populations
        return None

    def kinetics_matrix(self, l, m, F, wavelength_laser):
        # Compute the kinetics matrix expanded in l, m coefficients
        # Note that the order of the laser in F must be consistent with the
        # choices made in the fluorophore class.
        # In this case F[0] is the blue light.

        # l, m, F, and wavelength must be provided by outside.
        # l and m are the quantum numbers of SH, while F encodes the info about
        # photon flux and it's angular dependence.

        # Initialize arrays
        Feye = np.eye( (l.size) )
        K = np.zeros( (self.nspecies, self.nspecies, l.size, l.size))

        # Put the kinetic constants connecting the right species
        # K[1,0] is the kinetic constant for the proces 1 <- 0.
        nlasers = F.shape[0]
        nwavelengths = self.wavelength.size
        if self.nspecies == 4:
            self.fluorophore_type = 'rsFP_negative_4states'
            for i in np.arange(nlasers):
                for j in np.arange(nwavelengths):
                    if wavelength_laser[i] == self.wavelength[j]:
                        K[1,0] = K[1,0] + F[i] * self.cross_section_on[j]
                        K[3,2] = K[3,2] + F[i] * self.cross_section_off[j]
            K[0,1] = Feye / self.lifetime_on
            K[2,1] = Feye / self.lifetime_on  * self.quantum_yield_on_to_off
            K[2,3] = Feye / self.lifetime_off
            K[0,3] = Feye / self.lifetime_off * self.quantum_yield_off_to_on
        
        if self.nspecies == 6:
            self.fluorophore_type = 'rsFP_negative_6states'
            for i in np.arange(nlasers):
                for j in np.arange(nwavelengths):
                    if wavelength_laser[i] == self.wavelength[j]:
                        K[1,0] = K[1,0] + F[i] * self.cross_section_on[j]
                        K[4,3] = K[4,3] + F[i] * self.cross_section_off[j]
            K[0,1] = Feye / self.lifetime_on
            K[2,1] = Feye / self.lifetime_on  * self.quantum_yield_on_to_off
            K[3,2] = Feye / self.protonation_time_on
            K[3,4] = Feye / self.lifetime_off
            K[5,4] = Feye / self.lifetime_off * self.quantum_yield_off_to_on
            K[0,5] = Feye / self.deprotonation_time_off
        
        if self.nspecies == 8:
            self.fluorophore_type = 'rsFP_negative_8states'
            for i in np.arange(nlasers):
                for j in np.arange(nwavelengths):
                    if wavelength_laser[i] == self.wavelength[j]:
                        K[1,0] = K[1,0] + F[i] * self.cross_section_on[j]  # cis anionic
                        K[3,2] = K[3,2] + F[i] * self.cross_section_on[j]  # trans anionic
                        K[5,4] = K[5,4] + F[i] * self.cross_section_off[j]  # trans neutral
                        K[7,6] = K[7,6] + F[i] * self.cross_section_off[j]  # cis neutral
            # On-branch
            K[0,1] = Feye / self.lifetime_on
            K[2,1] = Feye / self.lifetime_on  * self.quantum_yield_on_to_off
            K[2,3] = Feye / self.lifetime_on
            K[0,3] = Feye / self.lifetime_on  * self.quantum_yield_trans_to_cis_anionic
            K[4,2] = Feye / self.protonation_time_on

            # Off-branch
            K[4,5] = Feye / self.lifetime_off
            K[6,5] = Feye / self.lifetime_off * self.quantum_yield_off_to_on
            K[6,7] = Feye / self.lifetime_off
            K[4,7] = Feye / self.lifetime_off * self.quantum_yield_cis_to_trans_neutral
            K[0,6] = Feye / self.deprotonation_time_off
        return K

    def starting_coeffs(self, l, m):
        # Compute the starting values for the population SH coefficients
        c0 = np.zeros( (self.nspecies, l.size) )
        c0[:,0] = self.starting_populations
        return c0

######################
# Illumination classes
######################
class SingleLaser:
    def __init__(self,
                 power_density,
                 polarization='x',
                 wavelength=488, 
                 numerical_aperture=1.4,
                 refractive_index=1.518):
        # Compute the photon flux from power density and wavelength
        self.power_density = power_density  # [W/cm2]
        self.wavelength = wavelength  # [nm]
        self.photon_flux = photon_flux(self.power_density,
                                       self.wavelength)  # [photons/cm2]

        # Polarization of the beam
        # It can be 'x', 'y', or 'c' for circular. 
        self.polarization = polarization

        # Numerical aperture of excitation light beam
        self.numerical_aperture = numerical_aperture

        # Refractive index of the immersion medium of the objective
        self.refractive_index = refractive_index

        # Time windows settings for compatibility
        self.modulation = np.array([[1]], dtype='float')
        self.time_windows = np.array([0], dtype='float')
        self.nwindows = 1
        self.time0 = 0
        return None
    
    def photon_flux_prod_coeffs(self, l, m, wigner_3j_prod_coeffs):
        # Compute the prod_matrix for the angular dependence of the photon flux.
        # When the photon flux is multiplied by the absorbtion cross-section
        # the kinetic rate of absortion is obtained.
        c_exc = np.zeros((1, l.size))
        c_exc[0] = na_corrected_linear_coeffs(l, m,
                polarization = self.polarization,
                numerical_aperture = self.numerical_aperture,
                refractive_index = self.refractive_index,
                ) * self.photon_flux
        
        # Return the matrix F as a 3 dimensional array, this will be helpfull
        # when more than one laser is used to excite the sample.
        F = np.zeros((1, l.size, l.size))
        F[0] = kinetic_prod_block(c_exc[0], wigner_3j_prod_coeffs)
        return F, c_exc

class ModulatedLasers:
    def __init__(self,
                 power_density,
                 polarization=['x', 'xy'],
                 wavelength=[405, 488],
                 modulation=[[1, 0], [1, 1]],
                 time_windows=[250e-9, 1e-3],
                 time0=0,
                 numerical_aperture=1.4,
                 refractive_index=1.518):
        # Compute the photon flux from power density and wavelength
        self.power_density = np.array(power_density)  # [W/cm2]
        self.wavelength = np.array(wavelength)  # [nm]
        self.photon_flux = photon_flux(self.power_density,
                                       self.wavelength)  # [photons/cm2]

        # Polarization of the beam
        # It can be 'x', 'y', or 'c' for circular. 
        self.polarization = polarization

        # Time Modulation properties
        self.modulation = np.array(modulation, dtype='float')
        self.time_windows = np.array(time_windows, dtype='float')
        self.nwindows = self.time_windows.shape[0]
        self.time0 = time0

        # Numerical aperture of excitation light beam
        # Refractive index of the immersion medium of the objective
        self.numerical_aperture = numerical_aperture
        self.refractive_index = refractive_index

        return None
    
    def photon_flux_prod_coeffs(self, l, m, wigner_3j_prod_coeffs):
        # Compute the prod_matrix for the angular dependence of the photon flux.
        # When the photon flux is multiplied by the absorbtion cross-section
        # the kinetic rate of absortion is obtained.

        # Initialize arrays for SH photoselection coefficients in c, and 
        # product coefficients in F.
        nlasers = np.size(self.polarization)
        c_exc = np.zeros((nlasers, l.size))
        F = np.zeros((nlasers, l.size, l.size))
        for i in np.arange(nlasers):
            c_exc[i] = na_corrected_linear_coeffs(l, m,
                    polarization = self.polarization[i],
                    numerical_aperture = self.numerical_aperture,
                    refractive_index = self.refractive_index,
                    ) * self.photon_flux[i]
            F[i] = kinetic_prod_block(c_exc[i], 
                        wigner_3j_prod_coeffs)
        
        # F is a 3 dimensional array. The firt index is for the lasers, the last
        # indexes are for lm of SH.
        return F, c_exc

def photon_flux(power_density, wavelength):    
    # Compute the photon flux from power density and wavelength.
    # The wavelength is used to compute the photon energy
    # photon_energy = h*c/wavelength
    # The wavelength is coverted in meters.
    light_speed = 299792457  # [m/s]
    h_planck = 6.62607015e-34  # [J*s]
    photon_energy = h_planck * light_speed / (wavelength*1e-9)  # [J]
    photon_flux = power_density / photon_energy  # [photons/cm2]
    return photon_flux

# TODO Check that photon flux is conserved changing polarization state

def na_corrected_linear_coeffs(l, m,
                               polarization='x', 
                               numerical_aperture=1.4,
                               refractive_index=1.518,
                               norm='4pi'):
    # Compute Axelrod high-NA correction coefficients
    k = axelrod_correction_k(numerical_aperture, refractive_index)

    # Get the SH expansion coefficients for x, y, and z photoselections.
    cx = linear_light_matter_coeffs(l, m, polarization='x', norm=norm)
    cy = linear_light_matter_coeffs(l, m, polarization='y', norm=norm)
    cz = linear_light_matter_coeffs(l, m, polarization='z', norm=norm)

    # Compute the SH expansion coefficients
    # We are assuming that the propagation direction is z.
    # If z polarization is chosen, the propagation direction is assumed as x.
    if polarization == 'x':
        c = k[0]*cz + k[1]*cy + k[2]*cx

    if polarization == 'y':
        c = k[0]*cz + k[1]*cx + k[2]*cy

    if polarization == 'z':
        c = k[0]*cx + k[1]*cy + k[2]*cz

    if polarization == 'xy':
        c = k[0]*cz + (k[1]+k[2])/2*cy + (k[1]+k[2])/2*cx

    # c are the SH coefficients of the expanded angular function.
    return c

def axelrod_correction_k(numerical_aperture, refractive_index):
    # Normalized Axelrod correction cartesian coefficients from [Fisz2005] eq.(2).
    # For a linear polarized beam:
    # - k[0] refers to the propagation direction
    # - k[1] refers to the perpendicular direction
    # - k[2] refers to the parallel direction
    # In the same paper the correction for the spherical coordinates problem are
    # discussed. It will be useful when adding the photoswithing angle.
    # Note that k[0] + k[1] + k[2] = 1. In other words, for a randomly oriented
    # sample in linear excitation, the total number of excited molecules is the
    # same regardless of the NA.

    # The maximum angle of the focalized beam with respect to the
    # propagation direction of the beam is computed.
    # NA = n*sin(theta), thus theta = arcsin(NA/n).
    max_ray_angle = np.arcsin(numerical_aperture/
                              refractive_index)

    k = np.zeros(3)
    costh = np.cos(max_ray_angle)
    k[0] = 1/6  * (2 -3*costh             +costh**3) / (1 - costh)
    k[1] = 1/24 * (1 -3*costh +3*costh**2 -costh**3) / (1 - costh)
    k[2] = 1/8  * (5 -3*costh   -costh**2 -costh**3) / (1 - costh)
    return k

def linear_light_matter_coeffs(l, m, polarization='x', 
                                    norm='4pi'):
    # This function computes the spherical harmonics expansion for linear 
    # light matter interactions photoselection functions.
    # For example, the function (r dot x)^2 is expanded and has only three non 
    # zero coefficients, in particular l=0, m=0 and l=2, m=0,2. r and x are 
    # versors. In this function few polarizations are expressed in terms of
    # their SH coefficients.
    c = np.zeros(l.shape)
    if polarization == 'x':
        c[np.logical_and(l==0, m==0)] = 1/3
        c[np.logical_and(l==2, m==0)] = -np.sqrt(1/45)
        c[np.logical_and(l==2, m==2)] = np.sqrt(1/15)
    if polarization == 'y':
        c[np.logical_and(l==0, m==0)] = 1/3
        c[np.logical_and(l==2, m==0)] = -np.sqrt(1/45)
        c[np.logical_and(l==2, m==2)] = -np.sqrt(1/15)
    if polarization == 'z':
        c[np.logical_and(l==0, m==0)] = 1/3
        c[np.logical_and(l==2, m==0)] = np.sqrt(4/45)
    if polarization == 'xy':  # circular polarization in the xy plane
        c[np.logical_and(l==0, m==0)] = 1/3
        c[np.logical_and(l==2, m==0)] = -np.sqrt(1/45)
    if norm == 'ortho':
        c = c/np.sqrt(4*np.pi)

    # Normalize the coefficients so the total integral of the function is 1
    # c = c/c[0]  # For now it is not on, because we loose the probability
    # connection
    return c

###################
# Detection classes
###################
class PolarizedDetection:
    def __init__(self,
                 polarization=['x', 'y'],
                 numerical_aperture=1.4,
                 refractive_index=1.518,
                 ):
        # Numerical aperture of the collection objective
        self.numerical_aperture = numerical_aperture

        # Refractive index of the immersion medium
        self.refractive_index = refractive_index

        # Polarization directions of the detection channels
        self.polarization = polarization
        return None

    def detector_coeffs(self, l, m):
        # Compute the SH expansion for the angular collection functions of the
        # detectors
        ndetectors = np.size(self.polarization)
        c = np.zeros((ndetectors, l.size))

        for i in np.arange(ndetectors):
            c[i] = na_corrected_linear_coeffs(l, m,
                polarization=self.polarization[i],
                numerical_aperture=self.numerical_aperture,
                refractive_index=self.refractive_index,
                )
        return c


################################################################################
# Engine of the rotational-diffusion and kinetics solver
################################################################################

def quantum_numbers(lmax):
    # Generate arrays with quantum numbers l and m
    l = np.array([])
    m = np.array([])
    for li in np.arange(lmax+1):
        # This is the way pyshtools is exporting vectors of sh expansions
        l = np.concatenate((l, np.ones(2*li+1)*li))
        m = np.concatenate((m, np.arange(li+1)))
        m = np.concatenate((m, np.arange(-1,-li-1,-1)))
    l = np.int32(l)
    m = np.int32(m)
    return l, m

def clebsch_gordan_prod_3darray(l, m):
    # 3D array with all the clebsh gordan coefficients for sh multiplication
    n = l.size
    cgp = np.zeros([n, n, n])
    for i in np.arange(n):
        for j in np.arange(n):
            for k in np.arange(n):
                    cgp[i,j,k] = np.sqrt( (2*l[i] + 1) * (2*l[j] + 1) / 
                                        ( np.pi*4    * (2*l[k] + 1) ) ) * (
                                sf.clebsch_gordan(l[i], 0,
                                                l[j], 0,
                                                l[k], 0) * 
                                sf.clebsch_gordan(l[i], m[i],
                                                l[j], m[j],
                                                l[k], m[k])
                                )
    
    # Constant due to normalizatino issues
    cgp = cgp*np.sqrt(4*np.pi) 
    return cgp

def wigner_3j_prod_3darray(l, m):
    # 3D array with all the coefficient for sh multiplication
    # Here we use Wigner3j symbols from sht, it is 10-20 time faster than
    # clebsch_gordan_prod_3darray(l, m).
    
    # Preliminary computations
    n = l.size
    lmax = np.uint(np.max(l))
    w3jp = np.zeros([n, n, n])

    # Limit calculation to allowed l1 indexes
    # This optimization uses the fact that l1 can assume only limite values 
    # in light-matter interaction.
    # For linear interaction only l1=0 and l1=2 are allowed.
    l1_allowed = np.logical_or(l==0, l==2)
    l1_allowed_indexes = np.arange(n)[l1_allowed]

    for i in l1_allowed_indexes:
        for j in np.arange(n):
            # Get all the quantum numbers for 1 and 2
            l1 = l[i]
            m1 = m[i]
            l2 = l[j]
            m2 = m[j]

            # Compute quantum numbers for 3
            m3 = -m1-m2 # requirement for w3j to be non zero

            # Compute Wigner3j symbols
            # w3j0 = w3jcalc.calculate(l1, l2, 0, 0) 
            # w3j1 = w3jcalc.calculate(l1, l2, m1, m2)
            w3j1 =  wigner_3j_all_l(l, l1, l2, m3, m1, m2, lmax)
            w3j0 =  wigner_3j_all_l_m0(l, l1, l2, lmax)

            # Compute SH product coefficients
            w3jp[i,j,:] = np.sqrt( (2*l1 + 1) * (2*l2 + 1) * (2*l + 1) / 
                                    (np.pi*4) ) * w3j0 * w3j1 *(-1)**np.float(m3) 
            #NOTE: (-1)** factor is necessary to match the result obtained
            # with clebsh-gordan coefficients. I am not sure why it is the case.
    
    # Constant due to normalizatino issues
    # Normalization of SH: https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html
    w3jp = w3jp*np.sqrt(4*np.pi)
    return w3jp

def wigner_3j_all_l(l, l1, l2, m3, m1, m2, lmax):
    # Compute all Wigner 3j symbols for a set of l3 at l1,l2,m3,m1,m2

    # Compute the coefficients 
    # https://shtools.github.io/SHTOOLS/pywigner3j.html
    w3j, l3min, l3max = sht.utils.Wigner3j(l1, l2, m3, m1, m2)
    l3 = np.arange(l3min, l3max+1)
    l3 = l3[l3<=lmax]  # Restrict to values below lmax 
    w3j = w3j[np.arange(l3.size)]

    # Index of l, m vector of SHTools
    # https://shtools.github.io/SHTOOLS/pyshcilmtovector.html
    i = np.uint(1.5-np.sign(-m3)/2)  # The sign minus is necessary, m3 = -M
    k = np.uint(l3**2+(i-1)*l3+np.abs(m3))

    # Final array of l.size with Wigner 3j coefficients
    w3jl = np.zeros(l.size)
    w3jl[k] = w3j
    return w3jl

def wigner_3j_all_l_m0(l, l1, l2, lmax):
    # Compute all Wigner 3j symbols for a set of l3 at l1,l2,m3,m1,m2

    # Compute the coefficients 
    # https://shtools.github.io/SHTOOLS/pywigner3j.html
    w3j, l3min, l3max = sht.utils.Wigner3j(l1, l2, 0, 0, 0)
    l3 = np.arange(l3min, l3max+1)
    l3 = l3[l3<=lmax]  # Restrict to values below lmax 
    w3j = w3j[np.arange(l3.size)]

    # Create array of size l.size adding w3j symbols for all l values
    w3jl = np.zeros(lmax+1)
    w3jl[l3] = w3j
    l3 = np.arange(lmax+1)
    w3jl = w3jl[l]
    return w3jl

def kinetic_prod_block(kvec, cgp):
    # Create a kinetic constant block for multiplication using cg coeffs.
    kblock = np.transpose(cgp, [2, 1, 0]).dot(kvec)
    return kblock

def isotropic_diffusion_block(l, m, diffusion_coefficient):
    # Make the diffusion diagonal block for isotropic rotational diffusion
    ddiag = -l*(l+1)*diffusion_coefficient
    dblock = np.zeros([l.size, l.size])
    np.fill_diagonal(dblock, ddiag)
    return dblock

def isotropic_diffusion_matrix(l, m, diffusion_coefficient, nspecies):
    # Make all the diffusion matrices for all the species assuming isotropic
    # rotational diffusion and the same rotational diffusion coefficient for
    # every species.
    D = np.zeros( (nspecies, l.size, l.size) )
    D[0] = isotropic_diffusion_block(l, m, diffusion_coefficient)
    for i in np.arange(1,nspecies):
        D[i] = D[0]
    return D

def diffusion_kinetics_matrix(D, K):
    # Create the full kinetic diffusion matrix expanded in sh
    nspecies = D.shape[0]
    M = np.zeros(K.shape)

    for i in range(nspecies):
        for j in range(nspecies):
            # Add the rotational diffusion blocks on the diagonal
            if i==j:
                M[i,i] = D[i]
            # Add the kinetics blocks in all the slots
            M[i,j] = M[i,j] + K[i,j]

    # Add on the diagonal the kinetics contributions that deplete the states
    for i in range(nspecies):
        for j in range(nspecies):
            if i != j:
                M[i,i] = M[i,i] - M[j,i]
    return M

def kinetics_diffusion_matrix_lmax(Dvec, Kmatrix, lmax):
    # Create the full kinetic diffusion matrix expanded in sh staring 
    # from scratch.
    # In this routine the expansion in l,m are included.
    assert len(Dvec) == len(Kmatrix)

    # Compute quantum number arrays and cg coefficients
    # The main simplification is that linear light matter interaction 
    # limits l1 to 0 and 2.
    l, m = quantum_numbers(lmax)
    #TODO Compute only cg coefficients that are needed, partially done
    # cgp = clebsch_gordan_prod_3darray(l, m)
    cgp = wigner_3j_prod_3darray(l, m)

    ncoeff = l.size  # number of spherical harmonics expansion coefficients
    nspecies = len(Dvec)  # number of species

    M = np.array( [ [np.zeros([ncoeff,ncoeff])]*nspecies ]*nspecies )
    for i in range(nspecies):
        for j in range(nspecies):
            if i == j:
                M[i,j] = isotropic_diffusion_block(l, m, Dvec[i])
            # Create blocks of the kinetic matrix rates
            else:
                k = Kmatrix[i][j]
                if np.size(k) == 1:
                    M[i,j] = np.eye(ncoeff).dot(k)
                else:
                    M[i,j] = kinetic_prod_block(k, cgp)

    for i in range(nspecies):
        for j in range(nspecies):
            if i != j:
                M[i,i] = M[i,i]-M[j,i]
    return M

def solve_evolution(M, c0, time):
    # Analitically solve the diffusion-kinetic problem by matrix exp of M
    nspecies = c0.shape[0]
    ncoeffs = c0.shape[1]
    c = np.zeros((nspecies, ncoeffs, time.size))

    # Variables with unique index for species and spherical harmonics
    c0 = c0.flatten()
    M = np.transpose(M, axes=[0,2,1,3])
    M = np.reshape(M, [nspecies*ncoeffs, nspecies*ncoeffs])
    c = np.reshape(c, [nspecies*ncoeffs, time.size])

    #TODO: Optimize matrix multiplication for only c0 values different from zero
    # Coefficients of starting conditions that are zeros
    # Useful for optimizing matrix multiplication
    # Most of the coefficient are zero because of simmetry and could be removed
    # from the problem.

    # Diagonalize M and invert eigenvector matrix
    # This will optimize the computation of matrix exponentiation.
    L, U = np.linalg.eig(M)
    Uinv = np.linalg.inv(U)

    # The following is equivalent to:
    #
    # for i in np.arange(time.size):
    #     ci = np.matmul(U, np.diag(np.exp(L*time[i])))
    #     ci = np.matmul(ci, Uinv)
    #     ci = ci.dot(c0)
    #     c[:,i] = ci
    #
    # Doing the matrix multiplication transposed reduces all the costly full
    # matrix multiplications to only matrix times vectors.
    # Depending on matrix size this version is several times faster.

    ci0 = c0.dot(Uinv.T)
    UTr = U.T  # Transpose the matrix only once
    for i in np.arange(time.size):
        ci = np.matmul(ci0, np.diag(np.exp(L*time[i])))
        ci = np.matmul(ci, UTr)
        c[:,i] = np.real(ci)  # Discard small rounding error complex parts

    # Reshape the coefficients into a 3D array, separating the coefficients
    # for each species.
    c = np.reshape(c, (nspecies, ncoeffs, time.size))
    return c, L, U

def remove_odd_coeffs(M, c0, l):
    l_sel = l%2 == 0
    c0 = c0[:, l_sel]
    M = M[:, :, :, l_sel, l_sel]
    return M, c0

def add_odd_coeffs_zeros(M, c0, c, l):
    c0_out = np.zeros((c0.shape[0], l.size))
    c_out = np.zeros((c0.shape[0], l.size, c0.shape[2]))
    M_out = np.zeros((M.shape[0], M.shape[1], M.shape[2], l.size, l.size))
    l_sel = l%2 == 0
    c0_out[:,l_sel] = c0
    c_out[:,l_sel,:] = c
    M_out[:,:,:,l_sel,l_sel] = M
    return M_out, c0_out, c_out

################################################################################
# Auxiliary functions
################################################################################

def make_angles(lmax):
    # Create angles for computing your functions of theta and phy
    cos2 = sht.SHGrid.from_zeros(lmax=lmax, kind='real')
    # cos2 = sht.shclasses.DHComplexGrid.from_zeros(lmax=2)
    theta = cos2.lats()
    phi = cos2.lons()
    omega = np.meshgrid(theta,phi)
    omega[0] = omega[0].transpose()*np.pi/180
    omega[1] = omega[1].transpose()*np.pi/180
    return omega

def make_grid(farray, lmax):
    # Make a sht grid from array data and expand in sh up to lmax
    fgrid = sht.SHGrid.from_zeros(lmax=lmax, kind='real')
    # fgrid = sht.SHGrid.from_array(farray)
    fgrid.data = farray
    fcilm = fgrid.expand(normalization='4pi')
    fvec = sht.shio.SHCilmToVector(fcilm.coeffs, lmax)

    # Remove small coefficients from errors
    fvec[np.abs(fvec)<1e-6] = 0 
    return fgrid, fvec, fcilm

def vecs2grids(cvecs):
    # Convert an array of vectors of SH coefficients to an array of SH grids
    cgridarray = []
    ntimes = cvecs.shape[1]
    for i in range(ntimes):
        cvec = cvecs[:,i]
        cgrid = vec2grid(cvec)
        cgridarray.append(cgrid)
    return cgridarray

def vec2grid(cvec, lmax=32):
    # Convert a vector of SH coefficients to an SH grid
    cilm = sht.shio.SHVectorToCilm(cvec)
    carr = sht.expand.MakeGridDH(cilm, sampling=2, extend=1, norm=1, lmax=lmax) # 4 orthonorm, 1 4pi
    cgrid = sht.shclasses.SHGrid.from_array(carr)
    return cgrid


################################################################################
# Plotting functions
################################################################################

def plot_prob_sphere(grid):
    prob = grid.data
    omega = make_angles(grid.lmax)
    cmap = plt.cm.Blues

    x = np.cos(omega[0]) * np.cos(omega[1])
    y = np.cos(omega[0]) * np.sin(omega[1])
    z = np.sin(omega[0])

    ax = plt.subplot(111, projection='3d')
    ax._axis3don = False
    ax.plot_surface(x, y, z,
                    cstride=1,
                    rstride=1,
                    facecolors=cmap(prob))
    ax.view_init(elev=30, azim=150)
    ax.set_xticks(np.linspace(-1,1,5))
    ax.set_yticks(np.linspace(-1,1,5))
    ax.set_zticks(np.linspace(-1,1,5))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    return

def plot_proj(grid, clims=[0, 1], cmap=plt.cm.Blues):
    prob = grid.data

    if clims == []:
        clims = [np.min(prob), np.max(prob)]

    im1 = plt.imshow(prob, extent=(0, 2, -0.5, 0.5),
                    vmin=clims[0], vmax=clims[1],
                    cmap=cmap, alpha=1)
    plt.xlabel('$\phi$ $(\pi)$')
    plt.ylabel('$\\theta$ $(\pi)$')
    plt.xticks([0, 0.5, 1, 1.5, 2])
    plt.yticks([0.5, 0.25, 0, -0.25, -0.5])
    plt.colorbar(ax=plt.gca())
    return


if __name__ == "__main__":
    from codetiming import Timer

    # rsEGFP2 = NegativeSwitcher(cross_section_on_blue=1e-10,
    #                            lifetime_on=3.6e-9,
    #                            quantum_yield_on_to_off=0.001,
    #                            diffusion_coefficient=1/(6*100e-6) )

    # tau = 100e-6 # us anisotropy decay
    # D = 1/(tau*6) # Hz
    # yield_off = 0.001 
    # tau_off = 80e-6 # time constant off switching
    # tau_on_exc = 3.6e-9 # lifetime of excited state
    # k21 = 1 / (tau_off * yield_off)
    # k12 = 1 / tau_on_exc
    # k32 = (1/tau_off) / yield_off

    # lmax = 8
    # omega = make_angles(lmax)
    # k21a = (np.sin(omega[0])**2) * k21
    # k21grid, k21c, k21cilm = make_grid(k21a, lmax)
    # plot_proj(k21grid, clims=[])

    # # Test product using cg coeff
    # l, m = quantum_numbers(lmax)
    # t = Timer()
    # t.start()
    # cgp = clebsch_gordan_prod_3darray(l, m)
    # t.stop()
    # # k21prod = kinetic_prod_block(k21c, cgp)
    # # kp = np.cos(omega[0])**2 * np.cos(omega[1]+np.pi)**2
    # # kpgrid, kpc = make_grid(kp, lmax)
    # # k21kpc = k21prod.dot(kpc)
    # # k21kpcilm = sht.shio.SHVectorToCilm(k21kpc)
    # # k21kparray = sht.expand.MakeGridDH(k21kpcilm, sampling=2)
    # # k21kpgrid = sht.shclasses.SHGrid.from_array(k21kparray)
    # # k21kpgrid.plot3d()

    # # Array with all the diffusion tensors/scalar for every specie.
    # # In this case every specie diffuse with the same rate.
    # Dvec = [D, D, D]

    # # Array with kinetic constants connecting the states.
    # Kmatrix = [[   0, k12, 0],
    #            [k21a,   0, 0],
    #            [   0, k32, 0]]


    # # # Simplified kinetic scheme
    # # # Array with all the diffusion tensors/scalar for every specie.
    # # # In this case every specie diffuse with the same rate.
    # # Dvec = [D, D]

    # # # Array with kinetic constants connecting the states.
    # # Kmatrix = [[   0,   0],
    # #            [k21a,   0]]

    # t = Timer()
    # t.start()
    # M = kinetics_diffusion_matrix(Dvec, Kmatrix, lmax)
    # t.stop()

    # # initial conditions
    # c0a = omega[0]*0 +1
    # c0grid, c0vec, c0cilm = make_grid(c0a, lmax)
    # c0 = np.zeros(((lmax+1)**2,3))
    # c0[:,0] = c0vec
    # # p0_1 = 1 + np.cos(omega[0]) * 0
    # # p0_1grid, c0_1 = make_grid(p0_1, lmax)
    # # c0[:,0] =c0_1

    # time = np.logspace(-11,-3,128)
    # # time = np.linspace(0,1e-3,1000)
    # t.start()
    # c, L, U =solve_evolution(M, c0, time)
    # t.stop()
    # # plt.imshow(np.real(Dblock))
    
    # Uinv = np.linalg.inv(U)
    # ci = np.matmul(U, np.diag(np.exp(L*100e-6)))
    # ci = np.matmul(ci, Uinv)
    # plt.imshow(np.real(ci))


    # t.start()
    # np.linalg.inv(U)
    # t.stop()

    # cplot = c[:,0,:]
    # cplotgrid = vecs2grids(cplot)

    # plt.figure()
    # plt.plot(time, cplot.T)
    # plt.xscale('log')

    # plot_proj(cplotgrid[100])

    # # cplotgrid[10].plot()
    # # plt.imshow(M)
    # # plt.show()
    # # plt.imshow(cplotgrid[20].data, vmin=0, vmax=1)