# ebvpy
A python function to calculate the **E(B-V)** reddening value for a source at a specific location on the sky using the Schlegel et al. (1998) dust maps. 

This function module provides a function for extracting the E(B-V) extinction/reddening value from the 
Schlegel et al. 1998, ApJ, 500, 525: http://adsabs.harvard.edu/abs/1998ApJ...500..525S. 

The maps are split into SGP and NGP and in the paper a prescription is given for converting between galactic longititude (l) and latitude (b) coordinates and x and y pixel coordinates of the maps. This conversion is described as:
<a href="http://www.codecogs.com/eqnedit.php?latex=x&space;=&space;2048&space;\times&space;\sqrt{1&space;-&space;n*\sin&space;(b)}&space;\cos&space;(l)&space;&plus;&space;2047.5&space;\\&space;y&space;=&space;-2048&space;\times&space;n&space;\sqrt{1&space;-&space;n*\sin&space;(b)}&space;\sin&space;(l)&space;&plus;&space;2047.5" target="_blank"><img src="http://latex.codecogs.com/gif.latex?x&space;=&space;2048&space;\times&space;\sqrt{1&space;-&space;n*\sin&space;(b)}&space;\cos&space;(l)&space;&plus;&space;2047.5&space;\\&space;y&space;=&space;-2048&space;\times&space;n&space;\sqrt{1&space;-&space;n*\sin&space;(b)}&space;\sin&space;(l)&space;&plus;&space;2047.5" title="x = 2048 \times \sqrt{1 - n*\sin (b)} \cos (l) + 2047.5 \\ y = -2048 \times n \sqrt{1 - n*\sin (b)} \sin (l) + 2047.5" /></a>

where n=+1 for the north galactic pole map (i.e. b > 0) and n = -1 for the south galactic pole map (i.e. b < 0).

There is also a function for calculating the magnitude correction, A_b in a given band from the A_b/E(B-V) values provided in Table 9 of Schlegel et al. (1998) or a user defined value. 

The corrected magnitude is then therefore m_corr,b = m_obs,b - A_b. 

An example use of the functions is as follows:

  `>>> from ebvpy import *`
  
  `>>> ebv = calc_ebv(ra=[160.23, 150.24, 26.346], dec=[-32.457, -4.578, 67.32])`
  
  `>>> print ebv`
  
  `   array([ 0.07547796,  0.0263974 ,  1.47947371])`
  
  `>>> calc_color_correction('u', ebv)`
  
  `   array([ 0.31995106,  0.11189858,  6.27148906])`
  


