# Quantum Dot Functions

import numpy as np

def quantum_yield():
	return 0.90+0.9*np.random.random()

def generatePDF(x,mean,variance):
	numerator = np.exp(-(x-mean)**2/(2*variance))
	denominator = np.sqrt(2*np.pi*variance)
	return numerator / denominator

def absorptionSpectra(x,mean,variance):

	return generatePDF(x,mean,variance)

def emissionSpectra(x,mean,variance):

	return generatePDF(x,mean,variance)

def emissionAngles():

	xi = np.random.normal()
	gamma_i = 2*np.pi*xi
	chi = 2*xi-1
	theta_e = np.arccos(np.sign(chi)-chi)

	return np.array(gamma_e,theta_e)
