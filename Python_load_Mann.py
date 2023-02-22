# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:10:51 2023
This is an example of how to load in a Mann box file and plot the contour and corresponding PSD. 
Use with care! I make many errors!
@author: cgrinde
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

# Load routines stolen from Wind Energy Toolbox on Gitlab: 
#    https://gitlab.windenergy.dtu.dk/toolbox/WindEnergyToolbox/-/blob/master/wetb/wind/turbulence/mann_turbulence.py

def load(filename, N=(32, 32)):
    """Load mann turbulence box

    Parameters
    ----------
    filename : str
        Filename of turbulence box
    N : tuple, (ny,nz) or (nx,ny,nz)
        Number of grid points

    Returns
    -------
    turbulence_box : nd_array

    Examples
    --------
    >>> u = load('turb_u.dat')
    """
    data = np.fromfile(filename, np.dtype('<f'), -1)
    if len(N) == 2:
        ny, nz = N
        nx = len(data) / (ny * nz)
        assert nx == int(nx), "Size of turbulence box (%d) does not match ny x nz (%d), nx=%.2f" % (
            len(data), ny * nz, nx)
        nx = int(nx)
    else:
        nx, ny, nz = N
        assert len(data) == nx * ny * \
            nz, "Size of turbulence box (%d) does not match nx x ny x nz (%d)" % (len(data), nx * ny * nz)
    return data.reshape(nx, ny * nz)


#Example! Make sure to change numbers according to your setup!
#Defining size of box
n1=4096
n2=32
n3=32

Lx=6142.5
Ly=180
Lz=180

umean=9

deltay=Ly/(n2-1)
deltax=Lx/(n1-1)
deltaz=Lz/(n3-1)
deltat=deltax/umean

time=np.arange(deltat, n1*deltat+deltat, deltat)

# Load in the files and reshape them into 3D
u=load('sim1.bin',  N=(n1, n2, n3))

ushp = np.reshape(u, (n1, n2, n3)) 


# Plot a countour
fig,ax=plt.subplots(1,1)
cp = ax.contourf(ushp[1000,:,:])
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
#ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.show()


# Picking a point on the velocity plane
sig=ushp[:,20,20]


fs=1/(time[1]-time[0])

#Compute and plot the power spectral density.
# Check out this site for inputs to the Welch
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html 
f, Pxx_den = signal.welch(sig, fs, nperseg=1024)
fig,ax=plt.subplots(1,1)
plt.loglog(f, Pxx_den)
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD ')
plt.show()