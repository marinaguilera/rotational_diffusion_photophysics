# %%
import pyshtools as sht
import spherical as sf
import numpy as np
import starss_kinetics_difusion_solution as kd

# %%

w3j, lmin, lmax = sht.utils.Wigner3j(3,3,2,1,0)

print(w3j)
print(lmin)
print(lmax)

# %%
lmax = 10
j2_max = lmax
j3_max = lmax
calc3j = sf.Wigner3jCalculator(j2_max, j3_max)

# Presumably, the following is inside some loop over j2, j3, m2, m3  
j2 = 5
j3 = 2
m2 = 0
m3 = 0
w3j = calc3j.calculate(j2, j3, m2, m3)  
m1 = - m2 - m3  
for j1 in range(max(abs(j2-j3), abs(m1)), j2+j3+1):  
    w3j[j1]  # This is the value of the 3-j symbol written above 


# %%
#%%
lmax=5
calc3j = sf.Wigner3jCalculator(lmax, lmax)

l1 = 1
l2 = 2
m1 = 1
m2 = -0

# Compute 3 quantum numbers
m3 = -m1-m2 # requirement for w3j to be non zero
i = np.int((1-np.sign(m3))/2)
l3 = np.arange(max(abs(l1-l2), abs(m3)), l1+l2+1)

w3j0 = calc3j.calculate(l1, l2, 0, 0) 
print(w3j0)
w3j1 = calc3j.calculate(l1, l2, m1, m2)
print(w3j1)
print(l3)

w3j1, l3min1, l3max1 = sht.utils.Wigner3j (l1, l2, m3, m1, m2)
w3j0, l3min0, l3max0 = sht.utils.Wigner3j (l1, l2, 0, 0, 0)


# %%
l, m = kd.quantum_numbers(lmax)
w3j = kd.wigner_3j_all_l(l,l1,l2,m3,m1,m2)
w3j

# %%
l3 = 2
sf.clebsch_gordan(l1, m1, l2, m2, l3, -m3)
# %%
