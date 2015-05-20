# ebvpy
A python function to calculate the **E(B-V)** reddening value for a source at a specific location on the sky using the Schlegel et al. (1998) dust maps. 

Once the E(B-V) value is derived, to calculate the change to the magnitude in a given band that must be made you must look up an appropriate A/E(B-V) value in the literature. 

A_v = m_obs - m_int 

A_v = E(B-V) * A/E(B-V)

