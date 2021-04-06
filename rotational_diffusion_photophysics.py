
import numpy as np
import matplotlib.pyplot as plt
import spherical as sf 
import pyshtools as sht 

# This module solves analytically the diffusion and kinetics
# in STARSS experiments.

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
                                    (np.pi*4) ) * w3j0 * w3j1 *(-1)**(2-np.mod(m3,2)) 
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

def kinetic_prod_block(kvec, cgp):
    # Create a kinetic constant block for multiplication using cg coeffs.
    n = kvec.size
    kblock = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            # Get slice of CG coefficients for the total l and m
            cgi = cgp[:,j,i].transpose()
            kblock[i,j] = cgi.dot(kvec)
    return kblock

def diffusion_block(D, l, m):
    # Make the diffusion diagonal block for isotropic rotational diffusion
    ddiag = -l*(l+1)*D
    dblock = np.zeros([l.size, l.size])
    np.fill_diagonal(dblock, ddiag)
    return dblock

def kinetics_diffusion_matrix(Dvec, Kmatrix, lmax):
    # Create the full kinetic diffusion matrix expanded in sh
    #TODO: this function doesn't work
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

    #TODO: Create M as a 2d list of matrices, this will help with the computation
    # of the full matrix.

    M = np.array( [ [np.zeros([ncoeff,ncoeff])]*nspecies ]*nspecies )
    for i in range(nspecies):
        for j in range(nspecies):
            if i == j:
                M[i,j] = diffusion_block(Dvec[i], l, m)
            # Create blocks of the kinetic matrix rates
            else:
                k = Kmatrix[i][j]
                if np.size(k) == 1:
                    M[i,j] = np.eye(ncoeff).dot(k)
                else:
                    kgrid, kvec, kcilm = make_grid(k, lmax)
                    M[i,j] = kinetic_prod_block(kvec, cgp)

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
        c[:,i] = ci

    # Reshape the coefficients into a 3D array, separating the coefficients
    # for each species.
    c = np.reshape(c, (nspecies, ncoeffs, time.size))
    return c, L, U

# Auxiliary functions
def vecs2grids(cvecs):
    # Convert an array of vectors of SH coefficients to an array of SH grids
    cgridarray = []
    ntimes = cvecs.shape[1]
    for i in range(ntimes):
        cvec = cvecs[:,i]
        cgrid = vec2grid(cvec)
        cgridarray.append(cgrid)
    return cgridarray

def vec2grid(cvec):
    # Convert a vector of SH coefficients to an SH grid
    cilm = sht.shio.SHVectorToCilm(cvec)
    carr = sht.expand.MakeGridDH(cilm, sampling=2, extend=1, norm=1, lmax=32) # 4 orthonorm, 1 4pi
    cgrid = sht.shclasses.SHGrid.from_array(carr)
    return cgrid

# Plotting functions
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

    fig, ax1 = plt.subplots()
    im1 = plt.imshow(prob, extent=(0, 2, -0.5, 0.5),
                    vmin=clims[0], vmax=clims[1],
                    cmap=cmap, alpha=1)
    plt.xlabel('$\phi$ $(\pi)$')
    plt.ylabel('$\\theta$ $(\pi)$')
    plt.xticks([0, 0.5, 1, 1.5, 2])
    plt.yticks([0.5, 0.25, 0, -0.25, -0.5])
    plt.colorbar(ax=ax1)
    return


if __name__ == "__main__":
    from codetiming import Timer

    tau = 100e-6 # us anisotropy decay
    D = 1/(tau*6) # Hz
    yield_off = 0.001 
    tau_off = 80e-6 # time constant off switching
    tau_on_exc = 3.6e-9 # lifetime of excited state
    k21 = 1 / (tau_off * yield_off)
    k12 = 1 / tau_on_exc
    k32 = (1/tau_off) / yield_off

    lmax = 8
    omega = make_angles(lmax)
    k21a = (np.sin(omega[0])**2) * k21
    k21grid, k21c, k21cilm = make_grid(k21a, lmax)
    plot_proj(k21grid, clims=[])

    # Test product using cg coeff
    l, m = quantum_numbers(lmax)
    t = Timer()
    t.start()
    cgp = clebsch_gordan_prod_3darray(l, m)
    t.stop()
    # k21prod = kinetic_prod_block(k21c, cgp)
    # kp = np.cos(omega[0])**2 * np.cos(omega[1]+np.pi)**2
    # kpgrid, kpc = make_grid(kp, lmax)
    # k21kpc = k21prod.dot(kpc)
    # k21kpcilm = sht.shio.SHVectorToCilm(k21kpc)
    # k21kparray = sht.expand.MakeGridDH(k21kpcilm, sampling=2)
    # k21kpgrid = sht.shclasses.SHGrid.from_array(k21kparray)
    # k21kpgrid.plot3d()

    # Array with all the diffusion tensors/scalar for every specie.
    # In this case every specie diffuse with the same rate.
    Dvec = [D, D, D]

    # Array with kinetic constants connecting the states.
    Kmatrix = [[   0, k12, 0],
               [k21a,   0, 0],
               [   0, k32, 0]]


    # # Simplified kinetic scheme
    # # Array with all the diffusion tensors/scalar for every specie.
    # # In this case every specie diffuse with the same rate.
    # Dvec = [D, D]

    # # Array with kinetic constants connecting the states.
    # Kmatrix = [[   0,   0],
    #            [k21a,   0]]

    t = Timer()
    t.start()
    M = kinetics_diffusion_matrix(Dvec, Kmatrix, lmax)
    t.stop()

    # initial conditions
    c0a = omega[0]*0 +1
    c0grid, c0vec, c0cilm = make_grid(c0a, lmax)
    c0 = np.zeros(((lmax+1)**2,3))
    c0[:,0] = c0vec
    # p0_1 = 1 + np.cos(omega[0]) * 0
    # p0_1grid, c0_1 = make_grid(p0_1, lmax)
    # c0[:,0] =c0_1

    time = np.logspace(-11,-3,128)
    # time = np.linspace(0,1e-3,1000)
    t.start()
    c, L, U =solve_evolution(M, c0, time)
    t.stop()
    # plt.imshow(np.real(Dblock))
    
    Uinv = np.linalg.inv(U)
    ci = np.matmul(U, np.diag(np.exp(L*100e-6)))
    ci = np.matmul(ci, Uinv)
    plt.imshow(np.real(ci))


    t.start()
    np.linalg.inv(U)
    t.stop()

    cplot = c[:,0,:]
    cplotgrid = vecs2grids(cplot)

    plt.figure()
    plt.plot(time, cplot.T)
    plt.xscale('log')

    plot_proj(cplotgrid[100])

    # cplotgrid[10].plot()
    # plt.imshow(M)
    # plt.show()
    # plt.imshow(cplotgrid[20].data, vmin=0, vmax=1)
