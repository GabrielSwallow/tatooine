#!/usr/bin/env python3

import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt

n=246

data = pluto.Pluto('import/011122_1/')
sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
print(np.argwhere(np.isnan(sigma)))
#sigma[np.isnan(sigma)] = 1.0
R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
sig_r = data.user_def_parameters['SIGMA_REF']

plt.figure()
plt.plot(R[1,:-1], np.mean(sigma, axis=0)/ sig_r  )
plt.plot(R, R**-0.5, '-')
#plt.ylim(0,0.1)
print(np.max(np.mean(sigma, axis=0)/ data.units['density'] ))
plt.xlim(0,20)
plt.xlabel(r'radius [a_bin]')
plt.ylabel(r'$\Sigma \, [g/cm^2]$')
plt.title('Kep 47')
plt.show()

# plt.figure()
# v_phi = data.primitive_variable('vx2', n)[0,:,:]
# plt.plot(R[1,:-1], np.mean(v_phi, axis=0))
# #plt.ylim(0,1)
# plt.xlim(0,20)
# plt.xlabel(r'radius [a_bin]')
# plt.ylabel(r'$v_{\phi} [T_{bin}]$')
# plt.title('')
# plt.show()

#%%

# for n in range(0,15):
#     print(n)
#     sigma = data.primitive_variable('rho', n)[0,:,:] #* data.units['density']
#     print(np.argwhere(np.isnan(sigma)))
#     #sigma[np.isnan(sigma)] = 1.0
#     R, Phi = np.meshgrid(data.grid['faces'][0], data.grid['faces'][1])
#     sig_r = data.user_def_parameters['SIGMA_REF']

#     plt.figure()
#     plt.plot(R[1,:-1], np.mean(sigma, axis=0)/ sig_r  )
#     plt.plot(R, R**-0.5, '-')
#     #plt.ylim(0,0.1)
#     print(np.max(np.mean(sigma, axis=0)/ data.units['density'] ))
#     plt.xlim(0,20)
#     plt.xlabel(r'radius [a_bin]')
#     plt.ylabel(r'$\Sigma \, [g/cm^2]$')
#     plt.title('Kep 47')
#     plt.show()

