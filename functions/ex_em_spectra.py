import numpy as np

def ex_spectrum(alpha, beta, peak):
    return lambda x: 2 ** (-((np.log(-alpha * (-(1/alpha) + x - peak))) ** 2) / (np.log(beta) ** 2))

def em_spectrum(alpha, beta, peak):
    return lambda x: 2 ** (-((np.log(-alpha * (-(1/alpha) - x + peak))) ** 2) / (np.log(beta) ** 2))

def find_ext_coeff(ex_spectrum, max_epsilon):
    return lambda x: (np.log(10) * max_epsilon) * ex_spectrum(x)
