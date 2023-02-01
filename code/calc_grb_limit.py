""" Calculate the GRB upper limit """

import numpy as np
import vals

# limits are the same
fluence = 1E-6
flux = 3E-7

Eiso = fluence * (4*np.pi*vals.dL_cm**2)
Liso = flux * (4*np.pi*vals.dL_cm**2)
print(Eiso, Liso)

