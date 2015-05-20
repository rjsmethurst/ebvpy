"""
This function module provides a function for extracting the E(B-V) extinction/reddening value from the 
Schlegel et al. 1998, ApJ, 500, 525: http://adsabs.harvard.edu/abs/1998ApJ...500..525S. 

The maps are split into SGP and NGP and in the paper a prescription is given for converting between galactic longititude (l) 
and latitude (b) coordinates and x and y pixel coordinates of the maps. This conversion is described as:

X=sqrt(1-n*sin(b))*cos(l)*2048 + 2047.5

Y=-n*sqrt(1-n*sin(b))*sin(l)*2048 + 2047.5

where n=+1 for the north galactic pole map (i.e. b > 0) and n = -1 for the south galactic pole map (i.e. b < 0). 

There is also a function for calculating the magnitude correction, A_b in a given band from the A_b/E(B-V) values provided 
in Table 9 of Schlegel et al. (1998) or a user defined value. 

The corrected magnitude is then therefore m_corr,b = m_obs,b - A_b. 

"""

import numpy as np 
from matplotlib import pyplot as plt 
import os.path
import sys

from astropy.io import fits 
from astropy import units as u
from astropy.coordinates import SkyCoord

import wget

def calc_xy_pos(l,b):
	""" Given galactic coordinate positions (l,b) i.e. longititude and latitude, this function calculates the appropriate [x,y] pixel positions 
	of the Schlegel et al. (1998) dust maps. 

	INPUTS
	=-=-=-
	:l:
		galactic longitude, in units of radians. Numpy array of dimensions (N,)
	:b:
		galactic latitude, in units of radians. Numpy array of dimensions (N,)

	OUTPUTS
	=-=-=-=

	:x:
		column pixel index position of the Schlegel et al. (1998) dust map. Numpy array with dimensions (N,) of the input array
	:y:
		row pixel index position of the Schlegel et al. (1998) dust map. Numpy array with dimensions (N,) of the input array

	"""
	n = np.array(b > 0).astype(int)
	n[n==0] = -1
	x = 2048 * ((1-(n*np.sin(b)))**(0.5))*np.cos(l) + 2047.5
	y = (-1)*2048 * n * ((1-(n*np.sin(b)))**(0.5))*np.sin(l) + 2047.5
	return np.round(x, 0).astype(int), np.round(y,0).astype(int)


def calc_ebv(ra, dec):
    """ Given sky coordinates (ra, dec) this function calculates the corresponding E(B-V) value on that part of the sky using the Schlegel et al. (1998) dust maps.

	INPUTS
	=-=-=-
	:ra:
		right ascension of sources, in units of degrees. Numpy array of dimensions (N,). Cannot take singular value, must be in list/array.
	:dec:
		declination of sources, in units of degrees. Numpy array of dimensions (N,). Cannot take singular value, must be in list/array.
	OUTPUTS
	=-=-=-=

	:ebv:
		E(B-V) reddening value of source coordinates provided. Numpy array with dimensions (N,) of the input array
	

	"""
    if os.path.exists('SFD_dust_4096_sgp.fits') == False:
		wget.download('http://nebel.rc.fas.harvard.edu/mjuric/lsd-data/sfd-dust-maps/SFD_dust_4096_sgp.fits', out='SFD_dust_4096_sgp.fits')
    else:
        pass
    if os.path.exists('SFD_dust_4096_ngp.fits') == False:
		wget.download('http://nebel.rc.fas.harvard.edu/mjuric/lsd-data/sfd-dust-maps/SFD_dust_4096_ngp.fits', out='SFD_dust_4096_ngp.fits')
    else:
		pass
    ngp = fits.open('SFD_dust_4096_ngp.fits')[0].data
    sgp = fits.open('SFD_dust_4096_sgp.fits')[0].data
    c = SkyCoord(ra=ra, dec=dec,  unit='deg')
    x, y = calc_xy_pos(c.galactic.l.radian, c.galactic.b.radian)
    ebv = np.zeros(len(ra))
    for row in range(len(ebv)):
        if c.galactic.b.radian[row] > 0:
            ebv[row] = ngp[y[row], x[row]]
        else:
            ebv[row] = sgp[y[row], x[row]]
    return ebv

def calc_color_correction(band, ebv):
    """ Given the E(B-V) value and a photometric bandpass this function calculates the corresponding magnitude correction using pre-defined or user defined A_b/E(B-V) values.

	INPUTS
	=-=-=-
	:band:
		str/dict
		Can be one of the five SDSS bands provided: "u",  "g", "r", "i" or "z" or 
		the user can specify their own band with A_b/E(B-V) value with a python dicionary format: e.g. {"nuv": 8.18}
	:ebv:
		E(B-V) reddening value of a single or multiple sources. Numpy array with dimensions (N,) of the input array

	OUTPUTS
	=-=-=-=

	:ab
		A_b magnitude correction for given bandpass, b. Numpy array with dimensions (N,) of the input array. 
	"""
    a_ebv = {"u": 4.239, "g":3.303, "r":2.285, "i":1.698, "z":1.263}
    if type(band) == dict:
        a_ebv.update(band)
        ab = a_ebv[b.keys()[0]]*ebv
    else:
        try:
            ab = a_ebv[band]*ebv
        except KeyError:
            print "Not a recognised band, try submitting your own band with value as a dictionary to the function e.g. {'nuv':8.18}. Currently available are "+str(a_ebv.keys())
            sys.exit()
    return ab

