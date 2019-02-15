import numpy as np
global n_vec
n_vec = np.array([0,0,1])

def reflect(i_vec,n_vec):
	# Requires np input vectors
	return i_vec - 2*(np.dot(i_vec,n_vec))*n_vec

def snells(n1,n2,theta_i):
	return np.arcsin(n1/n2*np.sin(theta_i))

def refract(n1,n2,i_vec,n_vec,theta_i):
	sin2_theta_t = (n1/n2)**2(1-np.cos(theta_i))
	return n1/n2*i_vec+(np.cos(theta_i)-np.sqrt(1-sin2_theta_t))*n_vec

def critAngle(n1,n2):
	if n1 > n2:
		return np.arcsin(n2/n1)
	else:
		return False

def rPerp(n1,n2,theta_i):
	theta_t = snells(n1,n2,theta_i)
	numerator = n1*np.cos(theta_i)-n2*np.cos(theta_t)
	denominator = n1*np.cos(theta_i)+n2*np.cos(theta_t)
	return (numerator/denominator)**2

def rPara(n1,n2,theta_i):
	theta_t = snells(n1,n2,theta_i)
	numerator = n2*np.cos(theta_i)-n1*np.cos(theta_t)
	denominator = n2*np.cos(theta_i)+n1*np.cos(theta_t)
	return (numerator/denominator)**2

def reflRatio(n1,n2,theta_i):
	if not critAngle(n1,n2):
		return 0.5*(rPerp(n1,n2,theta_i)+rPara(n1,n2,theta_i))
	else:
		return 1

def transRatio(n1,n2,theta_i):
	return 1 - reflRatio(n1,n2,theta_i)
