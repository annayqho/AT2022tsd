""" Synchrotron calculations 

Following Readhead synchrotron notes
"""

import numpy as np

# Constants
A = 1.586E12

def get_V(R):
    V = (4/3) * np.pi * R**3
    return V


def get_g(alpha,nu1,nu2):
    g_first_term = (2*alpha+2)/(2*alpha+1)
    g_second_term_top = nu2**(alpha+1/2)-nu1**(alpha+1/2)
    g_second_term_bottom = nu2**(alpha+1)-nu1**(alpha+1)
    g = g_first_term * (g_second_term_top / g_second_term_bottom)
    return g


def Beq(R, L, alpha, nu1, nu2):
    """ Equipartition magnetic field """
    V = get_V(R)
    g = get_g(alpha, nu1, nu2)
    B = (8*np.pi*A*g*L/V)**(2/7)
    return B


def Ueq(R, L, alpha, nu1, nu2):
    V = get_V(R)
    g = get_g(alpha, nu1, nu2)
    B = Beq(R, L, alpha, nu1, nu2)
    U = 2 * V*B**2/(8*np.pi) # + A*g*L*B**(-3/2)
    return U


def calc_gamma_e(B, obs_nu):
    # CONSTANTS
    c = 3E10
    q_e = 4.8E-10
    m_e = 9.1E-28

    nu_g = q_e * B / (2 * np.pi * m_e * c)
    gamma_e = np.sqrt((obs_nu/nu_g))
    return gamma_e


if __name__=="__main__":
    # Lorentz factor
    Gamma = 1

    # Input: our values
    R = 9E11 * Gamma**2
    L = 1E43
    alpha = -1.6
    nu1 = 1E13
    nu2 = 1E15
    obs_nu = 1E15 # Hz

    B = Beq(R, L, alpha, nu1, nu2)
    U = Ueq(R, L, alpha, nu1, nu2)
    gamma_e = calc_gamma_e(B, obs_nu)

    tcool = 6*np.pi*9E-28*3E10 / (6.65E-25 * B**2 * gamma_e)
    #print(tcool) # cooling time in hours
