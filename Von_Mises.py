import numpy as np

# Von mises function, takes as input:
# - z position
# - y position
# - shear flow in zy plane
# - shear center loaction
# - moments about y and z
# Stress in x direction is the sum of the stress caused by bending around both axis
# Ixx and Iyy must be known beforehand
# outputs Von mises stress
def Von_Mises(z, y, tau_yz, z_hat, My, Mz ):

    sigma_xx_z = My * (z - z_hat) / Iyy
    sigma_xx_y = Mz * y / Izz

    sigma_xx_total = sigma_xx_y + sigma_xx_z
    sigma_vm =  np.square(0.5 * ((sigma_xx_total**2) + ((-sigma_xx_total)**2)) + 3*(tau_yz ** 2))

    return sigma_vm

